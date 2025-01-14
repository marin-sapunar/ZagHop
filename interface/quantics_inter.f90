!--------------------------------------------------------------------------------------------------
! PROGRAM: AnalyticModelInterface
!> @author Cris Sanz Sanz, Graham Worth
!> @author Marin Sapunar, Ruđer Boškovrć Institute
!> @date June, 2017
!
! Interface from Zagreb SH code to Quantics operator library
!--------------------------------------------------------------------------------------------------
module shzagreb_inter
      use decimal
      use global
      use versions
      use constants
      
      use dvrdatmod
      use griddatmod
      use operdef
      use rddvrmod
      use rdopermod
      use iorst, only: rstinfo
      use dirdyn, only: ndoftsh,dercpdim,ndofddpes,ltshtrans,ndofdd,&
                  tshtransb,dbnrec,nactdim,natmtsh,tshxcoo0,ldbsave,&
                  lupdhes,lnactdb,lddrddb,ddtrajnum,num_gp
      use dirdyn, only: alloc_dirdyn,alloc_dddb,atnam
      use directdyn
      use potevalmod, only: calcdiab,calcdiabder
      use psidef, only: qcentdim,gwpdim,zcent,vdimgp,dimgp,ndimgp,zgp,nsgp,totgp,&
                        sbaspar,rsbaspar
      use openmpmod, only: lompqc
      
      use dd_db, only: dddb_gp,getdbnrec,preparedb
      use dbcootrans
      use channels
      use op2lib, only: subvxxdo1
      use xvlib, only: mvxxdd1, mvtxdd1


      implicit none

      real(dop), allocatable :: qcoo(:),qcoo1(:),xgp(:)
      real(dop), allocatable :: tempvec(:, :, :)
      real(dop), allocatable :: pesdia(:,:)
      real(dop), allocatable :: derdia(:,:,:)
      real(dop), allocatable :: derad(:,:,:)
      complex(dop), allocatable :: rotmatz(:,:)
      integer(long), allocatable :: point(:)
      real(dop), allocatable :: hops(:)
      integer :: gdof

contains

      subroutine shzagreb_run(step, xyz0, cstate, en, gra, nadvec)

      
      integer, intent(in) :: step
      real(dop), intent(in) :: xyz0(:, :)
      integer, intent(in) :: cstate
      real(dop), intent(out) :: gra(:, :)
      real(dop), intent(out) :: en(:)
      real(dop), intent(out) :: nadvec(:, :, :)
    

      integer :: i,j,ilbl,jlbl,chkdvr,chkgrd,chkpsi,chkprp,n,f,f1,m
      logical(kind=4) :: linwf
      real(dop), allocatable :: xyz(:, :)

      real(dop) :: time

      integer             :: ierr
      character(len=c5)   :: filename,string
      logical(kind=4)     :: check,lerr
      real(dop), external :: dlamch
      logical(kind=4), save :: initialized=.false.

      allocate(xyz, source=xyz0)
      open(ilog,file='quantics.log',status='unknown',position='append')

      if (.not. initialized) then
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
            if (lddrddb) then
               lupdhes = .true.
            else
               lupdhes = .false.
            endif

! if not using a DB, do not allocate large memory
            if (.not. (lddrddb)) then
               dbmemdim = 1  
               ldbsmall = .false.
            endif
         endif

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
            call alloc_dirdyn
         endif
         if (lddtrans) call alloc_dbcootrans

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
         close(ioper)

!-----------------------------------------------------------------------
! Prepare DB and allocate memory
!-----------------------------------------------------------------------
         if (ldddb) then
            call preparedb(1)
            call getdbnrec(dbnrec)
            call alloc_dddb
         endif

!-----------------------------------------------------------------------
! Read data needed by the operator
!-----------------------------------------------------------------------
         operfile=oname(1:olaenge)//'/oper'
         allocate(hops(hopsdim))
         chkdvr=2
         chkgrd=1
         call rdoper(hops,chkdvr,chkgrd)

!-----------------------------------------------------------------------
! DD needs to read restart file for info on how Shepard Interpolation is 
! being done
!-----------------------------------------------------------------------
         if (ldd) then
           filename=rname(1:rlaenge)//'/restart'
           ilbl=index(filename,' ')-1
           open(irst,file=filename,form='unformatted',status='old',iostat=ierr)
           if (ierr .ne. 0) then
              routine='SHzagreb_interface'
              ilbl=index(filename,' ')-1
              message = 'Cannot open file: '//filename(1:ilbl)
              call errormsg
           endif
           chkdvr=0
           chkgrd=0
           chkpsi=0
           chkprp=1
           call rstinfo(linwf,lerr,chkdvr,chkgrd,chkpsi,chkprp)
           close(irst)
         endif

!-----------------------------------------------------------------------
! qcentdim is needed in getddpes as dimension of Ndof (effectively 1GWP)
!-----------------------------------------------------------------------
         if (ldd) then
            gwpdim(1,1) = 1
            zcent(1,1) = 1
            gwpm(1)=.true.
            qcentdim = nspfdof(1)

!needed for getddpes (No. of configurations is 1)
            vdimgp(1,1) = 1
            dimgp(1,1) = 1
            ndimgp(1,1) = 1
            nsgp(1)=1
            totgp = 1
            zgp(1) = 1

! no. of NACTS
            if (nddstate .gt. 1) then
               lnactdb = .true.
               nactdim = nddstate*(nddstate-1)/2
            endif

! get trajectory number from directory name
            string=ddname
            call dirpath(string,ilbl)
            jlbl=ilbl-1
            ilbl=jlbl
            do 
               if (string(ilbl-1:ilbl-1) .eq. '.') exit
               ilbl=ilbl-1
            enddo
            read(string(ilbl:jlbl),*) ddtrajnum

         endif

! get no. of states and no. of dynamical coordinates
         nstate = gdim(feb)
         gdof = nspfdof(1)

! Allocate memory
         allocate(tempvec(gdof,nstate,nstate))
         allocate(qcoo(ndoftsh))
         allocate(qcoo1(maxdim))
         allocate(xgp(maxdim))
         allocate(pesdia(maxsta,maxsta))
         allocate(rotmatz(maxsta,maxsta))
         allocate(point(maxdim))
         allocate(derdia(maxsta,maxsta,maxdim))
         allocate(derad(maxsta,maxsta,maxdim))

         initialized = .true.
      endif ! Initialization done

      gra = 0.0_dop
      nadvec = 0.0_dop

! reform xyz -> qcoo (Quantics dynamical coordinates)
      if (ltshtrans) then
         call subvxxdo1(xyz,tshxcoo0,ndoftsh)
         call mvxxdd1(tshtransb,xyz,qcoo,maxdim,ndoftsh,gdof)
      else
         f=0
         do n=1,natmtsh
            do i=1,3
               f=f+1
               qcoo(f)=xyz(i,n)
            enddo
         enddo
      endif

! need to add frozen coordinates to qcoo 
! (it assumes coordinates are the centre of a GWP)
      f1 = 0
      qcoo1 = 0.0
      m=1
      do n=1,nspfdof(m)
         f=spfdof(n,m)
         qcoo1(f) = qcoo(n)
      enddo
! Add in any frozen coordinates
      do f=1,ndof
         if (basis(f) .eq. 19) qcoo1(f) = rpbaspar(1,f)
      enddo

! Initialise local DBs. Need to be in Cartesians.
      if (ldd .and. ldbsmall) then
         if (lddtrans) then
            call ddq2x(qcoo1,xgp)
         else
            xgp=qcoo1
         endif

         num_gp = 1
         call dddb_gp(dbnrec,xgp,num_gp)
      endif

! Perform calculation of QC.
      time=0.0d0
      if (ldd) call getddpes(time,qcoo,1,1)   

    !  if (dddiab .eq. 0) then
! Extract energies, gradients, nacts from adiabatic data
    !     call extradgra(cstate,en,gra,nadvec, &
    !          dbener(1:nddstate,1:nddstate,1), &
    !          dbgrad(1:ndofddpes,1:nddstate,1:nddstate,1), &
    !          dbdercp(1:ndofddpes,1:dercpdim,1))
    !  else
! PES matrix in adiabatic (en) and diabatic (pesdia) representations
! rotmatz is the ADT matrix (as a complex)
         call calcdiab(hops,en,pesdia,rotmatz,point,qcoo1,1)

! Matrix of gradients in adiabatic (derad) and diabatic (derdia).
         call calcdiabder(hops,derad,derdia,rotmatz,qcoo1,1)

! Convert forces to Cartesian and extract forces / nact
         call extrgra(cstate,en,gra,nadvec,derad)
   !   endif

      close(ilog)

      end subroutine shzagreb_run


!#######################################################################

      subroutine extrgra(sta,en,gra,nadvec,derad)

      implicit none

      integer(long), intent(in)  :: sta
      integer(long)              :: s,s1,f,f1,n,m
      real(dop), dimension(ndoftsh,nddstate,nddstate), intent(out) :: nadvec
      real(dop), dimension(ndof)                                   :: qnadvec
      real(dop), dimension(ndoftsh),intent(out)              :: gra
      real(dop), dimension(ndof)                             :: qgra
      real(dop), dimension(maxsta,maxsta,maxdim), intent(in) :: derad
      real(dop), dimension(maxsta), intent(in)               :: en
      real(dop) :: ediff


! extract gradient for present state and transform to Cartesian
! removing frozen coordinates
      qgra(:) = 0.0
      m = 1  ! only 1 mode in TSH
      do n=1,nspfdof(m)
         f=spfdof(n,m)
         qgra(n) = derad(sta,sta,f)
      enddo
    
      gra(:) = 0.0
      if (ltshtrans) then
         call mvtxdd1(tshtransb,qgra,gra,maxdim,nspfdof(1),ndoftsh)
      else
         do f=1,ndoftsh
            gra(f)=qgra(f)
         enddo
      endif

! extract nacts and transform to Cartesian
! set  for frozen coordinates to 0
      qnadvec(:) = 0.0_dop
      nadvec(:,:,:) = 0.0_dop
      m = 1  ! only 1 mode in TSH
      do s=1,nddstate
         do s1=s+1,nddstate
            do n=1,nspfdof(m)
               f=spfdof(n,m)
               qnadvec(n) = derad(s1,s,f)
            enddo

            if (ltshtrans) then
               call mvtxdd1(tshtransb,qnadvec,nadvec(1,s1,s),maxdim,&
                    nspfdof(1),ndoftsh)
            else
               do f=1,ndoftsh
                  nadvec(f,s1,s)=qnadvec(f)
               enddo
            endif

            ediff = en(s) - en(s1)
            if (abs(ediff) .lt. 1.0d-6) then
               if (ediff .lt. 0.0) then
                  ediff=-1.0d-6
               else
                  ediff=1.0d-6
               endif
            endif
            nadvec(:,s1,s) = nadvec(:,s1,s) / ediff

            nadvec(:,s,s1) = -nadvec(:,s1,s)
         enddo
      enddo
  
      end subroutine extrgra

!#######################################################################

      subroutine extradgra(sta,en,gra,nadvec,av,deriv1,dercp)

      implicit none

      integer(long), intent(in)  :: sta
      integer(long)              :: s,s1,f,f1,sdx
      real(dop), dimension(ndoftsh,nddstate,nddstate), intent(out) :: nadvec
      real(dop), dimension(ndof)                                   :: qnadvec
      real(dop), dimension(ndoftsh),intent(out)              :: gra
      real(dop), dimension(ndofddpes)                        :: tmpgra
      real(dop), dimension(ndofdd)                           :: qgra
      real(dop), dimension(maxsta), intent(out)              :: en
      real(dop), dimension(nddstate,nddstate), intent(in)       :: av
      real(dop), dimension(ndofddpes,nddstate,nddstate), intent(in)  :: deriv1
      real(dop), dimension(ndofddpes,dercpdim), intent(in)           :: dercp
      real(dop) :: ediff

      gra = 0.0
      nadvec = 0.0

      do s=1,nddstate
         en(s) = av(s,s)
      enddo

! extract gradient for present state and transform to Cartesian
! removing frozen coordinates
      f1 = 0
      do f=1,ndofddpes
         f1 = f1+1
         tmpgra(f1) = deriv1(f,sta,sta)
      enddo

! first transform from Cartesians (DB coordinates) to DD coordinates
      if (lddtrans) then
         call ddf2q(tmpgra,qgra)
      else
         qgra(1:ndofdd) = tmpgra(1:ndofdd)
      endif
    
! transform from DD coordinates to TSH coordinates
      if (ltshtrans) then
         call mvtxdd1(tshtransb,qgra,gra,maxdim,nspfdof(1),ndoftsh)
      else
         do f=1,ndoftsh
            gra(f)=qgra(f)
         enddo
      endif

! extract nacts and transform to Cartesian
! set  for frozen coordinates to 0
      nadvec = 0.0_dop
      do s=1,nddstate-1
         do s1=s+1,nddstate
            sdx = (s1-1)*(s1-2)/2+s

            f1 = 0
            do f=1,ndofddpes
               f1 = f1+1
               tmpgra(f1) = dercp(f,sdx)
            enddo

            if (lddtrans) then
               call ddf2q(tmpgra,qnadvec)
            else
               qnadvec(1:ndofdd) = tmpgra(ndofdd)
            endif

            if (ltshtrans) then
               call mvtxdd1(tshtransb,qnadvec,nadvec(1,s1,s),maxdim,&
                    nspfdof(1),ndoftsh)
            else
               do f=1,ndoftsh
                  nadvec(f,s1,s)=qnadvec(f)
               enddo
            endif

            ediff = en(s) - en(s1)
            if (abs(ediff) .lt. 1.0d-6) then
               if (ediff .lt. 0.0) then
                  ediff=-1.0d-6
               else
                  ediff=1.0d-6
               endif
            endif
            nadvec(:,s1,s) = nadvec(:,s1,s) / ediff

            nadvec(:,s,s1) = -nadvec(:,s1,s)
         enddo
      enddo
  
      end subroutine extradgra

end module shzagreb_inter

