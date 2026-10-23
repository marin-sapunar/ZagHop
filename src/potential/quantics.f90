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
! diabatic potential is reordered to the ZagHop convention after each evaluation.
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
        dddiab,ddupdate,lddsoc
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
        integer, allocatable :: state_map(:) !< Quantics index of each ZagHop state (n_state).
    contains
        procedure :: init => quantics_initialize
        procedure :: set_group_nstate => quantics_set_group_nstate
        procedure :: update_geometry => quantics_update_geometry
        procedure :: eval_diab => quantics_eval_diab
        !   procedure :: get_oscill => quantics_get_oscill
    end type quantics_interface


contains

    subroutine quantics_update_geometry(self, geometry, coordinates)
        class(quantics_interface), intent(inout) :: self
        real(dop), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates
        integer :: i, n, f

        if (coordinates /= self%classical_coordinates) then
            write(stderr, *) 'Error in quantics_mod, update_geometry subroutine.'
            write(stderr, *) '  Unsupported coordinates "'//coordinates//'". Expected "'//          &
            &                self%classical_coordinates//'".'
            stop
        end if

        if (size(geometry, 1) /= 1 .or. size(geometry, 2) /= self%ndof) then
            write(stderr, *) 'Error in quantics_mod, update_geometry subroutine.'
            write(stderr, *) '  Geometry should have dimensions (1, ndof).'
            stop
        end if

        self%q = geometry(1, :)

        ! Place dynamical coordinates at their Quantics DOF index (frozen DOFs are set in init).
        do n = 1, self%ndof
            self%gpoint(spfdof(n, 1)) = self%q(n)
        end do
        if (allocated(self%adiab_trans)) self%ldiab_trans = self%adiab_trans
    end subroutine quantics_update_geometry


    subroutine quantics_initialize(self)
        class(quantics_interface) :: self

        character(len=c5)   :: filename, string
        logical(kind=4) :: check
        logical(kind=4) :: lerr
        integer :: ilbl, ierr, icheck, f, i
        integer :: chkdvr, chkgrd
        !, string
        !integer :: ilbl, jlbl, ierr, chkdvr, chkgrd, chkpsi, chkprp
        

        !integer(long) :: check
        !logical(kind=4) :: lcheck, lerr
        real(dop), external :: dlamch

        self%classical_coordinates = 'quantics'

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

        self%n_state = nddstate
        self%ndof = nspfdof(1)
        allocate(self%diab_w(self%n_state, self%n_state))
        allocate(self%diab_dw(self%ndof, self%n_state, self%n_state))
        allocate(self%adiab_w(self%n_state, self%n_state))
        allocate(self%adiab_dw(self%ndof, self%n_state, self%n_state))
        allocate(self%group_adiab_w(self%n_state, self%n_state))
        allocate(self%group_adiab_dw(self%ndof, self%n_state, self%n_state))
        allocate(self%group_adiab_trans(self%n_state, self%n_state))
        allocate(self%q(self%ndof), source=0.0_dop)
        ! Identity mapping until the spin multiplicities are given in set_group_nstate.
        self%state_map = [(i, i = 1, self%n_state)]

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


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: quantics_set_group_nstate
    !> @brief Set state groups and the mapping from ZagHop to Quantics state ordering.
    !> @details
    !! The oper file contains no information on the multiplicities of the states, so they are
    !! taken from the input. The states in each spin block are assumed to be ordered as in Quantics,
    !! with all components of each state consecutive.
    !----------------------------------------------------------------------------------------------
    subroutine quantics_set_group_nstate(self, nstate, multiplicity)
        use evaluator_base_mod, only : base_set_group_nstate => set_group_nstate
        class(quantics_interface), intent(inout) :: self
        integer, intent(in) :: nstate(:) !< Number of states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.
        integer :: i, j, k, ci, offset

        if (sum(nstate * multiplicity) /= self%n_state) then
            write(stderr, *) 'Error in quantics_mod, set_group_nstate subroutine.'
            write(stderr, '(3x,a,i0,a)') 'Input nstate and multiplicity give ', &
            &                            sum(nstate * multiplicity), ' states.'
            write(stderr, '(3x,a,i0,a)') 'Quantics operator has ', self%n_state, ' states.'
            stop
        end if

        call base_set_group_nstate(self, nstate, multiplicity)

        ci = 1
        offset = 0
        do i = 1, size(nstate)
            do j = 1, multiplicity(i)
                do k = 1, nstate(i)
                    self%state_map(ci) = offset + (k - 1) * multiplicity(i) + j
                    ci = ci + 1
                end do
            end do
            offset = offset + nstate(i) * multiplicity(i)
        end do
    end subroutine quantics_set_group_nstate


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

        ! Reorder states to ZagHop convention and extract derivatives with respect to the
        ! dynamical coordinates.
        self%diab_w = pesdia(self%state_map, self%state_map)
        do n = 1, nspfdof(m)
            f = spfdof(n, m)
            self%diab_dw(n, :, :) = derdia(self%state_map, self%state_map, f)
        end do

        close(ilog)
    end subroutine quantics_eval_diab

!#######################################################################

    ! subroutine extrsoc(pesspdi,sovec, spinvec)
    !     implicit none

    !     integer(long)              :: s,s1,f,f1,n,m
    !     real(dop), dimension(:,:), intent(out) :: sovec
    !     real(dop), dimension(:,:), intent(in)  :: pesspdi
    !     integer, dimension(:), intent(in) :: spinvec

    !     if(size(spinvec,1).gt.0)then
    !         ! Keep only values between different multiplicities
    !         do s=1,size(spinvec,1)
    !             do s1=1,size(spinvec,1)
    !             if (spinvec(s) .ne. spinvec(s1)) then
    !                 sovec(s1,s) = pesspdi(s1,s)
    !             endif
    !             enddo
    !         enddo
    !     elseif(allocated(imultmap))then !Ask Graham, cause nddstate is the number of multiplicity blocks, no the multiplicity of all states
    !         ! Keep only values between different multiplicities
    !         do s=1,nddstate
    !             do s1=1,nddstate
    !             if (imultmap(s) .ne. imultmap(s1)) then
    !                 sovec(s1,s) = pesspdi(s1,s)
    !             endif
    !             enddo
    !         enddo
    !     endif

    ! end subroutine extrsoc


! !#######################################################################

!     subroutine extrgra(sta,en,gra,nadvec,derad)
!         integer(long), intent(in)  :: sta
!         integer(long)              :: s,s1,f,f1,n,m
!         real(dop), dimension(ndoftsh,nddstate,nddstate), intent(out) :: nadvec
!         real(dop), dimension(ndof)                                   :: qnadvec
!         real(dop), dimension(ndoftsh),intent(out)                    :: gra
!         real(dop), dimension(ndof)                                   :: qgra
!         real(dop), dimension(nddstate,nddstate,maxdim), intent(in)   :: derad
!         real(dop), dimension(nddstate), intent(in)                   :: en
!         real(dop) :: ediff


!         ! extract gradient for present state and transform to Cartesian
!         ! removing frozen coordinates
!         qgra(:) = 0.0
!         m = 1  ! only 1 mode in TSH
!         do n=1,nspfdof(m)
!             f=spfdof(n,m)
!             qgra(n) = derad(sta,sta,f)
!         enddo

!         gra(:) = 0.0
!         if (ltshtrans) then
!             call mvtxdd1(tshtransb,qgra,gra,maxdim,nspfdof(1),ndoftsh)
!         else
!             do f=1,ndoftsh
!                 gra(f)=qgra(f)
!             enddo
!         endif

!     ! extract nacts and transform to Cartesian
!     ! set  for frozen coordinates to 0
!         qnadvec(:) = 0.0_dop
!         nadvec(:,:,:) = 0.0_dop
!         m = 1  ! only 1 mode in TSH
!         do s=1,nddstate
!             do s1=s+1,nddstate
!                 if (imultmap(s) .ne. imultmap(s1)) cycle  ! ignore different spins

!                 do n=1,nspfdof(m)
!                 f=spfdof(n,m)
!                 qnadvec(n) = derad(s1,s,f)
!                 enddo

!                 if (ltshtrans) then
!                 call mvtxdd1(tshtransb,qnadvec,nadvec(1,s1,s),maxdim,&
!                     nspfdof(1),ndoftsh)
!                 else
!                 do f=1,ndoftsh
!                     nadvec(f,s1,s)=qnadvec(f)
!                 enddo
!                 endif

!                 ediff = en(s) - en(s1)
!                 if (abs(ediff) .lt. 1.0d-6) then
!                 if (ediff .lt. 0.0) then
!                     ediff=-1.0d-6
!                 else
!                     ediff=1.0d-6
!                 endif
!                 endif
!                 nadvec(:,s1,s) = nadvec(:,s1,s) / ediff

!                 nadvec(:,s,s1) = -nadvec(:,s1,s)
!             enddo
!         enddo

!     end subroutine extrgra

end module quantics_interface_mod

