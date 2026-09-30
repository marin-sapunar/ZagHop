!--------------------------------------------------------------------------------------------------
! MODULE: quantics_interface_mod
!> @author Cris Sanz Sanz, Graham Worth
!> @author Marin Sapunar, Ruđer Boškovrć Institute
!> @date June, 2017
!
! Interface from Zagreb SH code to Quantics operator library
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
        lupdhes,lnactdb,lddrddb,ddtrajnum,num_gp
    use dirdyn, only: alloc_dirdyn,alloc_dddb,atnam,nsmult,imultmap
    use directdyn
    use potevalmod, only: calcdiab,calcdiabder,calcvreps,calcpes,calcpesder
    use psidef, only: qcentdim,gwpdim,zcent,vdimgp,dimgp,ndimgp,zgp,nsgp,totgp,&
        sbaspar,rsbaspar
    use openmpmod, only: lompqc
    use lalib, only: simtranbd

    use dd_db, only: dddb_gp,getdbnrec,preparedb
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
        !   procedure :: get_oscill => quantics_get_oscill
    end type quantics_interface


contains

    subroutine quantics_update_geometry(self, geometry)
        class(quantics_interface), intent(inout) :: self
        real(dop), intent(in) :: geometry(:, :)
        integer :: i, n, f

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

        ! ! Initialise local DBs. Need to be in Cartesians.
        ! if (ldd .and. ldbsmall) then
        !     if (lddtrans) then
        !         call ddq2x(self%gpoint,xgp)
        !     else
        !         xgp=self%gpoint
        !     endif

        !     num_gp = 1
        !     call dddb_gp(dbnrec,xgp,num_gp)
        ! endif
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
            call alloc_dirdyn(ilog)
        endif

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
        ! Read data needed by the operator
        !-----------------------------------------------------------------------
        operfile=oname(1:olaenge)//'/oper'
        allocate(self%hops(hopsdim))
        chkdvr=2
        chkgrd=1
        call rdoper(self%hops,chkdvr,chkgrd)

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

        ! Frozen coordinates are fixed at the centre of their basis.
        allocate(self%gpoint(maxdim), source=0.0_dop)
        do f = 1, ndof
            if (basis(f) == 19) self%gpoint(f) = rpbaspar(1, f)
        end do

    end subroutine quantics_initialize


    subroutine quantics_eval_diab(self)
        class(quantics_interface), intent(inout) :: self
        integer :: n, f, s, s1
        integer, parameter :: m = 1
        integer :: izflag, nham
        real(dop) :: time
        integer(long), allocatable :: point(:)
        complex(dop), allocatable :: pesdiaz(:, :)
        real(dop), allocatable :: derdia(:, :, :)

        open(ilog,file='quantics.log',status='unknown',position='append')

        ! Perform calculation of QC.
        time=0.0d0
        if (ldd) call getddpes(time,self%q,1,1)

        ! Calculate only diabatic potential and its derivatives.
        nham=1
        allocate(point(maxdim), source=1)
        allocate(pesdiaz(nddstate, nddstate))
        allocate(derdia(nddstate, nddstate, maxdim))
        call calcpes(self%hops, self%diab_w, point, self%gpoint, pesdiaz, izflag, nham)
        if (izflag /= 0) then
            write(stderr, *) 'Error in quantics_mod, eval_diab subroutine.'
            write(stderr, *) '  Complex potentials are not supported.'
            stop
        end if
        call calcpesder(self%hops, derdia, self%gpoint, nham)

        ! Extract derivatives with respect to the dynamical coordinates.
        do s = 1, nddstate
            do s1 = 1, nddstate
                do n = 1, nspfdof(m)
                    f = spfdof(n, m)
                    self%diab_dw(n, s1, s) = derdia(s1, s, f)
                end do
            end do
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

