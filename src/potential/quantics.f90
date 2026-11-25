!--------------------------------------------------------------------------------------------------
! MODULE: quantics_interface_mod
!> @author Cris Sanz Sanz, Graham Worth
!> @author Marin Sapunar, Ruđer Boškovrć Institute
!> @date June, 2017
!
! Interface from Zagreb SH code to Quantics operator library
!
! States of higher multiplicity are ordered differently in Quantics and ZagHop. In Quantics, all
! components of each state are consecutive (T1(1), T1(2), T1(3), T2(1), ...), while ZagHop places
! each component in a separate group of states (T1(1), T2(1), ..., T1(2), T2(2), ...). The
! evaluator declares the Quantics order (component_order = 'by_state'), and the base class
! reorders the diabatic potential after each evaluation (set_full_hamiltonian).
!
! The classical coordinates are the dynamical DOFs of the operator (spfdof(:, 1)), as written by
! Quantics in the geom and veloc files. For DD in Cartesian coordinates, they are the Cartesian
! coordinates of the atoms (classical_coordinates = 'cartesian', geometry (3, natom)), and otherwise
! one value per DOF (classical_coordinates = 'quantics', geometry (1, ndof)). Quantics uses the same
! rule. The masses are the effective masses of the kinetic energy terms of the operator
! (tshdofmass). For DD in Cartesians, these are the atom masses (equal for the x, y, z of an atom).
!
! For direct dynamics (DD) operators, the potential is interpolated from a Quantics database (DB).
! The trajectory is treated as a single GWP by the Quantics DD routines. If the DB is not read-only,
! new quantum chemistry calculations are run and added to the DB as required at each step.
!--------------------------------------------------------------------------------------------------
module quantics_interface_mod
    use evaluator_diabatic_mod
  !  use decimal
    use global
    use versions
    use constants

    use dvrdatmod
    use griddatmod
    use operdef
    use rddvrmod
    use rdopermod
    use iorst, only: rstinfo
    use dirdyn, only: dercpdim,ndofddpes,&
        dbnrec,nactdim,natmtsh,ldbsave,&
        lupdhes,lnactdb,lddrddb,lddwrdb,ddtrajnum,num_gp,&
        dddiab,ddupdate,lddsoc,tshdofmass,ltshdofmass
    use dirdyn, only: alloc_dirdyn,alloc_dddb0,alloc_dddb,atnam,nsmult,imultmap
    use directdyn
    use potevalmod, only: calcdiab,calcdiabder,calcvreps,calcpes,calcpesder
    use psidef, only: qcentdim,gwpdim,zcent,vdimgp,dimgp,ndimgp,zgp,nsgp,totgp,&
        sbaspar,rsbaspar
    use openmpmod, only: lompqc
    use lalib, only: simtranbd

    use dd_db, only: dddb_gp,getdbnrec,preparedb,setdbparams
    use logintmod, only: logdd
    use dbcootrans
    use channels
    use op2lib, only: subvxxdo1
    use xvlib, only: mvxxdd1, mvtxdd1

    implicit none

    type, extends(diabatic_evaluator) :: quantics_interface
        real(dop), allocatable :: q(:) !< Dynamical coordinates (nspfdof(1)).
        real(dop), allocatable :: gpoint(:) !< All Quantics DOFs (maxdim), including frozen DOFs.
        logical :: initialized = .false.
        real(dop), allocatable :: hops(:) !< Operator data read from the oper file (hopsdim).
    contains
        procedure :: init => quantics_initialize
        procedure :: update_geometry => quantics_update_geometry
        procedure :: eval_diab => quantics_eval_diab
        procedure :: default_reference => quantics_default_reference
        procedure :: check_atoms => quantics_check_atoms
        procedure :: set_model_states => quantics_set_model_states
        procedure :: save_state => quantics_save_state
        procedure :: load_state => quantics_load_state
    end type quantics_interface


contains

    !> The origin of the dynamical coordinates of the operator, q = 0. There is no default
    !! reference for DD in Cartesian coordinates.
    subroutine quantics_default_reference(self, geom)
        class(quantics_interface), intent(in) :: self
        real(dop), allocatable, intent(out) :: geom(:, :)

        if (self%classical_coordinates == 'quantics') then
            allocate(geom(1, self%ndof), source=0.0_dop)
        end if
    end subroutine quantics_default_reference


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: quantics_check_atoms
    !> @brief Check the masses of the system against the masses of the operator.
    !> @details
    !! With one mass per DOF (ndim = 1), the masses are compared with the effective masses of the
    !! kinetic energy terms (tshdofmass). With one mass per atom (DD in Cartesians), they are
    !! compared with the mass of the x coordinate of each atom. Labels are not checked. Oper files written before
    !! tshdofmass was added are not checked. The tolerance allows for masses given with 6
    !! significant digits.
    !----------------------------------------------------------------------------------------------
    subroutine quantics_check_atoms(self, sym, mass)
        class(quantics_interface), intent(in) :: self
        character(len=*), intent(in) :: sym(:) !< Atom symbols.
        real(dop), intent(in) :: mass(:) !< Atom masses (amu).
        real(dop), parameter :: tol = 1.0e-5_dop
        real(dop), allocatable :: ref(:)
        integer :: i

        if (.not. ltshdofmass) return
        if (size(mass) == self%ndof) then
            ref = [(tshdofmass(spfdof(i, 1)), i = 1, self%ndof)]
        else if (self%classical_coordinates == 'cartesian' .and. size(mass) == self%ndof / 3) then
            ref = [(tshdofmass(spfdof(3 * i - 2, 1)), i = 1, self%ndof / 3)]
        else
            write(stderr, *) 'Error in quantics_mod, check_atoms subroutine.'
            write(stderr, '(3x,a,i0,a,i0,a)') 'The system has ', size(mass),                       &
            &   ' atoms. The operator has ', self%ndof, ' dynamical DOFs.'
            stop
        end if
        ref = ref / Da_me
        do i = 1, size(mass)
            if (abs(mass(i) - ref(i)) > tol * ref(i)) then
                write(stderr, *) 'Error in quantics_mod, check_atoms subroutine.'
                write(stderr, '(3x,a,i0,3a)') 'Mass of atom ', i, ' (', trim(sym(i)),              &
                &   ') does not match the operator.'
                write(stderr, '(3x,a,es18.10,a,es18.10,a)') 'Input: ', mass(i), ' amu, operator: ', &
                &   ref(i), ' amu.'
                stop
            end if
        end do
    end subroutine quantics_check_atoms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: quantics_set_model_states
    !> @brief Set the states of the model (see evaluator_set_model_states) and check that
    !! their number matches the number of states of the operator.
    !----------------------------------------------------------------------------------------------
    subroutine quantics_set_model_states(self, nstate, multiplicity, spin_states)
        class(quantics_interface), intent(inout) :: self
        integer, intent(in) :: nstate(:) !< Number of spin-free states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.
        character(len=*), intent(in) :: spin_states

        call evaluator_set_model_states(self, nstate, multiplicity, spin_states)
        if (self%model_states%n_state /= nddstate) then
            write(stderr, *) 'Error in quantics_mod, set_model_states subroutine.'
            write(stderr, '(3x,a,i0,a)') 'Input nstate, multiplicity and spin_states give ',     &
            &                            self%model_states%n_state, ' states.'
            write(stderr, '(3x,a,i0,a)') 'The operator has ', nddstate, ' states.'
            stop
        end if
    end subroutine quantics_set_model_states


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: quantics_update_geometry
    !> @brief Set the dynamical coordinates of the operator.
    !> @details
    !! Accepts one value per DOF ('quantics', shape (1, ndof)) and, for DD in Cartesians, the
    !! Cartesian coordinates of the atoms ('cartesian', shape (3, ndof/3)). The DOFs of each atom
    !! are x, y, z in the order of the atoms, so both give the same flat vector, and gradients and
    !! couplings need no conversion.
    !----------------------------------------------------------------------------------------------
    subroutine quantics_update_geometry(self, geometry, coordinates)
        class(quantics_interface), intent(inout) :: self
        real(dop), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates
        integer :: n

        select case(coordinates)
        case('quantics')
            if (size(geometry, 1) /= 1 .or. size(geometry, 2) /= self%ndof) then
                write(stderr, *) 'Error in quantics_mod, update_geometry subroutine.'
                write(stderr, *) '  Geometry should have dimensions (1, ndof).'
                stop
            end if
        case('cartesian')
            if (self%classical_coordinates /= 'cartesian') then
                write(stderr, *) 'Error in quantics_mod, update_geometry subroutine.'
                write(stderr, *) '  Cartesian coordinates are only available for DD in Cartesians.'
                stop
            end if
            if (size(geometry, 1) /= 3 .or. size(geometry, 2) /= self%ndof / 3) then
                write(stderr, *) 'Error in quantics_mod, update_geometry subroutine.'
                write(stderr, *) '  Geometry should have dimensions (3, ndof/3).'
                stop
            end if
        case default
            write(stderr, *) 'Error in quantics_mod, update_geometry subroutine.'
            write(stderr, *) '  Unsupported coordinates "'//coordinates//'". Expected "'//          &
            &                self%classical_coordinates//'".'
            stop
        end select

        self%q = reshape(geometry, [self%ndof])

        ! Place dynamical coordinates at their Quantics DOF index (frozen DOFs are set in init).
        do n = 1, self%ndof
            self%gpoint(spfdof(n, 1)) = self%q(n)
        end do
    end subroutine quantics_update_geometry


    subroutine quantics_initialize(self)
        class(quantics_interface) :: self

        character(len=c5)   :: filename, string
        logical(kind=4) :: check
        logical(kind=4) :: lerr
        integer :: ilbl, ierr, icheck, f
        integer :: chkdvr, chkgrd
        !, string
        !integer :: ilbl, jlbl, ierr, chkdvr, chkgrd, chkpsi, chkprp
        

        !integer(long) :: check
        !logical(kind=4) :: lcheck, lerr
        real(dop), external :: dlamch

        open(ilog,file='quantics.log',status='unknown',position='append')

        macheps = dlamch('P')

        string='../..'
        ilbl=5
        call abspath(string,ilbl)
        dname = string
        dlaenge = index(dname,' ')-1
        oname = string
        olaenge = index(oname,' ')-1
        rname = string
        rlaenge = index(rname,' ')-1

        ! turn off parallelisation of QC calcs (omp threads do separate trajs).
        lompqc=.false.

        !-----------------------------------------------------------------------
        ! get array dimensions
        !-----------------------------------------------------------------------
        inquire(irst,opened=check)
        if (check) close(irst)
        filename=rname(1:rlaenge)//'/restart'
        ilbl=index(filename,' ')-1
        open(irst,file=filename(1:ilbl),form='unformatted',status='old',&
        iostat=ierr)
        if (ierr .ne. 0) then
            routine='SHzagreb_interface'
            ilbl=index(filename,' ')-1
            message = 'Cannot open file: '//filename(1:ilbl)
            call errormsg
        endif
        call rdmemdim(irst)
        close(irst)

        !-----------------------------------------------------------------------
        ! For DD calculations, DB needs to be read into memory
        !-----------------------------------------------------------------------
        if (ldd) then
            ldbsave = .true.
            lupdhes = lddrddb
            ! If not using a DB, do not allocate large memory.
            if (.not. lddrddb) then
                dbmemdim = 1
                ldbsmall = .false.
            end if
        end if

        !-----------------------------------------------------------------------
        ! Allocate memory
        !-----------------------------------------------------------------------
        allocmemory=0
        call alloc_dvrdat
        call alloc_grddat
        call alloc_operdef
        if (ldd .or. ltraj) then
            allocate(gwpdim(1,1))
            allocate(zcent(1,1))
            allocate(vdimgp(1,1))
            allocate(dimgp(1,1))
            allocate(ndimgp(1,1))
            allocate(zgp(1))
            allocate(nsgp(1))
            allocate(rsbaspar(sbaspar,maxdim,1))
        endif
        ! Needed for DD and for operators using a DB.
        call alloc_dirdyn(ilog)

        !-----------------------------------------------------------------------
        ! Read system / DVR information
        !-----------------------------------------------------------------------
         filename=dname(1:dlaenge)//'/dvr'
         ilbl=index(filename,' ')-1
         open(idvr,file=filename,form='unformatted',status='old',iostat=ierr)
         if (ierr .ne. 0) then
            routine='SHzagreb_interface'
            ilbl=index(filename,' ')-1
            message = 'Cannot open file: '//filename(1:ilbl)
            call errormsg
         endif
         chkdvr=1
         call dvrinfo(lerr,chkdvr)
         close(idvr)

        !-----------------------------------------------------------------------
        ! Read data from oper file
        !-----------------------------------------------------------------------
         ddpath = ' '
         filename=oname(1:olaenge)//'/oper'
         ilbl=index(filename,' ')-1
         open(ioper,file=filename,form='unformatted',status='old',iostat=ierr)
         if (ierr .ne. 0) then
            routine='SHzagreb_interface'
            ilbl=index(filename,' ')-1
            message = 'Cannot open file: '//filename(1:ilbl)
            call errormsg
         endif
         chkdvr=1
         chkdvr=2
         chkgrd=1
         call operinfo(lerr,chkdvr,chkgrd)

        !----------------------------------------------------------------------- 
        ! read in coordinate transformation information
        !-----------------------------------------------------------------------
        if (lddtrans .or. ltshtrans) then
            call alloc_dbcootrans
            call rdddtrans(ioper)
        endif

        close(ioper)

        !-----------------------------------------------------------------------
        ! Open DB and find out what it contains before allocating memory
        !-----------------------------------------------------------------------
        if (ldddb) then
            call alloc_dddb0(ilog)
            if (ldd .or. lddmctdh) then
                call preparedb(zeroi)
            else
                call preparedb(twoi)
            endif
            call setdbparams
            call alloc_dddb(ilog)
        endif

        !-----------------------------------------------------------------------
        ! Read data needed by the operator (also loads the DB into memory)
        !-----------------------------------------------------------------------
        operfile=oname(1:olaenge)//'/oper'
        allocate(self%hops(hopsdim))
        chkdvr=2
        chkgrd=1
        call rdoper(self%hops,chkdvr,chkgrd)

        if (ldd) call init_direct_dynamics()

        close(ilog)
        self%initialized = .true.

        ! The multiplicities are not in the oper file. They are given to set_model_states from the
        ! input, and checked against the number of states of the operator.
        self%component_order = 'by_state'
        self%ndof = nspfdof(1)
        allocate(self%q(self%ndof), source=0.0_dop)

        ! Coordinates of the trajectory, chosen as in Quantics (prop_zagreb).
        if (ldd .and. ddcoo == 2) then
            self%classical_coordinates = 'cartesian'
        else
            self%classical_coordinates = 'quantics'
        end if

        ! Frozen coordinates are fixed at the centre of their basis.
        allocate(self%gpoint(maxdim), source=0.0_dop)
        do f = 1, ndof
            if (basis(f) == 19) self%gpoint(f) = rpbaspar(1, f)
        end do

    end subroutine quantics_initialize


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: init_direct_dynamics
    !> @brief Set up the Quantics DD routines for a single trajectory.
    !> @details
    !! The DD options used during the propagation (DB update and interpolation scheme) are read
    !! from the restart file. The trajectory is passed to getddpes as a single GWP.
    !----------------------------------------------------------------------------------------------
    subroutine init_direct_dynamics()
        character(len=c5) :: filename
        logical(kind=4) :: linwf, lerr
        integer :: ilbl, ierr, chkdvr, chkgrd, chkpsi, chkprp

        filename=rname(1:rlaenge)//'/restart'
        ilbl=index(filename,' ')-1
        open(irst,file=filename(1:ilbl),form='unformatted',status='old',iostat=ierr)
        if (ierr .ne. 0) then
            routine='SHzagreb_interface'
            message = 'Cannot open file: '//filename(1:ilbl)
            call errormsg
        endif
        chkdvr=0
        chkgrd=0
        chkpsi=0
        chkprp=1
        call rstinfo(linwf,lerr,chkdvr,chkgrd,chkpsi,chkprp)
        close(irst)
        call logdd

        if (any(gwpm(2:nmode))) then
            write(stderr, *) 'Error in quantics_mod, init_direct_dynamics subroutine.'
            write(stderr, *) '  Direct dynamics requires all dynamical coordinates in the first mode.'
            stop
        end if
        if (lsocbas .or. lddsoc) then
            write(stderr, *) 'Error in quantics_mod, init_direct_dynamics subroutine.'
            write(stderr, *) '  Spin-orbit coupling is not supported with direct dynamics.'
            stop
        end if
        if (dddiab == 0 .and. nddstate > 1) then
            write(stderr, *) 'Error in quantics_mod, init_direct_dynamics subroutine.'
            write(stderr, *) '  Direct dynamics requires a diabatised DB (dddiab /= 0).'
            stop
        end if
        if (ddupdate == 1) then
            write(stderr, *) 'Error in quantics_mod, init_direct_dynamics subroutine.'
            write(stderr, *) '  Time controlled DB updates (ddupdate = 1) are not supported.'
            stop
        end if

        ! The trajectory is a single GWP in a single configuration.
        gwpdim(1,1) = 1
        zcent(1,1) = 1
        gwpm(1) = .true.
        qcentdim = nspfdof(1)
        vdimgp(1,1) = 1
        dimgp(1,1) = 1
        ndimgp(1,1) = 1
        nsgp(1) = 1
        totgp = 1
        zgp(1) = 1

        if (nddstate > 1) then
            lnactdb = .true.
            nactdim = nddstate*(nddstate-1)/2
        end if
    end subroutine init_direct_dynamics


    subroutine quantics_eval_diab(self)
        class(quantics_interface), intent(inout) :: self
        integer :: n, f
        integer, parameter :: m = 1
        integer :: izflag, nham
        real(dop) :: time
        integer(long), allocatable :: point(:)
        real(dop), allocatable :: pesdia(:, :)
        complex(dop), allocatable :: pesdiaz(:, :)
        real(dop), allocatable :: derdia(:, :, :)
        real(dop), allocatable :: xgp(:)
        real(dop), allocatable :: dw(:, :, :)

        open(ilog,file='quantics.log',status='unknown',position='append')

        if (ldd) then
            ! Set up the local DB around the current geometry. Distances are in DB coordinates.
            if (ldbsmall) then
                allocate(xgp(ndofddpes))
                if (lddtrans) then
                    call ddq2x(self%gpoint, xgp)
                else
                    xgp = self%gpoint(1:ndofddpes)
                end if
                num_gp = 1
                call dddb_gp(dbnrec, xgp, num_gp)
            end if

            ! Perform new QC calculations if needed and add them to the DB. Time is only used for
            ! time controlled updates, which are rejected in init_direct_dynamics.
            time=0.0d0
            if (lddwrdb .or. .not. lddrddb) call getddpes(time,self%q,1,1)
        end if

        ! Calculate only diabatic potential and its derivatives.
        nham=1
        allocate(point(maxdim), source=1)
        allocate(pesdia(nddstate, nddstate))
        allocate(pesdiaz(nddstate, nddstate))
        allocate(derdia(nddstate, nddstate, maxdim))
        call calcpes(self%hops, pesdia, point, self%gpoint, pesdiaz, izflag, nham)
        if (izflag /= 0) then
            write(stderr, *) 'Error in quantics_mod, eval_diab subroutine.'
            write(stderr, *) '  Complex potentials are not supported.'
            stop
        end if
        call calcpesder(self%hops, derdia, self%gpoint, nham)

        ! Derivatives with respect to the dynamical coordinates. The base class reorders the states
        ! to the ZagHop convention.
        allocate(dw(self%ndof, nddstate, nddstate))
        do n = 1, nspfdof(m)
            f = spfdof(n, m)
            dw(n, :, :) = derdia(:, :, f)
        end do
        call self%set_full_hamiltonian(pesdia, dw)

        close(ilog)
    end subroutine quantics_eval_diab

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: quantics_save_state
    !> @brief Add the state of the evaluated point to a JSON object, with the coordinates of the
    !! operator in the object 'quantics'.
    !> @details
    !! A DD database is not saved. Points added to it during the run are in the database files.
    !----------------------------------------------------------------------------------------------
    subroutine quantics_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        class(quantics_interface), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: qobj

        call diabatic_save_state(self, json, obj)
        call json%create_object(qobj, 'quantics')
        call json%add(qobj, 'q', self%q)
        call json%add(qobj, 'gpoint', self%gpoint)
        call json%add(obj, qobj)
    end subroutine quantics_save_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: quantics_load_state
    !> @brief Restore the state of an evaluated point saved by save_state.
    !----------------------------------------------------------------------------------------------
    subroutine quantics_load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child
        class(quantics_interface), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the state.
        type(json_value), pointer :: qobj

        call diabatic_load_state(self, json, obj)
        qobj => json_get_child(json, obj, 'quantics')
        call json_get_required(json, qobj, 'q', self%q, self%ndof)
        call json_get_required(json, qobj, 'gpoint', self%gpoint, size(self%gpoint))
    end subroutine quantics_load_state

end module quantics_interface_mod
