module mod_gridYpara_01y2
   implicit none

   !---------------------------------------------------------------------------------------!
   ! Determination of  grid for the AB + C  in a grid on R=R2, r=R1,gam                    !
   !                                01 + 2
   !
   !  input_grid: determines grid
   !  paralelizacion: determines index matrices and distributes quantities
   !_______________________________________________________________________________________!
   save
   ! * constants
   complex*16, parameter :: zero = dcmplx(0.d0,0.d0)
   complex*16, parameter :: zeye = dcmplx(0.d0,1.d0)
   real*8,parameter :: pi = dacos(-1.d0)
   real*8,parameter :: conve1 = 1.197d-4
   real*8,parameter :: convl = 0.52917726D0
   real*8,parameter :: convm = 0.182288853D4
   real*8,parameter :: conve = 4.55633538D-6
   real*8,parameter :: eV2cm = 8065.5d0
   real*8,parameter :: hbr = convl/dsqrt(convm*conve/conve1)
   real*8,parameter :: au2eV = 27.21113957d0
   real*8,parameter :: zot2au = 0.0380646230576441d0 !(1.d0/conve1)*conve
   real*8,parameter :: zot2eV = 1.03578177070099d0
   real*8,parameter :: hartree2cm = 219474.625

   real*8,parameter :: cluz_au = 137.d0
   real*8,parameter :: epsilon0_au = 0.07957747d0
   real*8,parameter :: Aconstant_au = 1.d0/(3.d0*pi*(cluz_au**3)*epsilon0_au) ! = 1/(3 pi hbar^4 Epsilon_0 c^3)
   real*8,parameter :: CSconstant_au = 1.d0/(cluz_au*epsilon0_au)    ! = 1/(hbar^2 Epsilon_0 c)

   ! *     public ! for input data in namelist, grid, basis

   ! *     grid and basis data in namelist
   integer :: npun1, npun1min, npun2, nangu, nangu2
   integer :: Jtot, iparity, inc, nelecmax, iommin, iommax, j0
   integer :: nvini, nvmax, jini, jmax
   integer :: nvref, jref, iomref, ielecref, ncan
   integer :: nproc, idproc

   ! *     write options
   integer :: iwrt_pot, iwrt_wvp, iwrt_reac_distri, n1plot, n2plot, nangplot

   ! *     process

   integer :: iphoto, nelec_bnd
   real*8  :: photonorm

   ! *     products states analysis

   integer :: iprod, nviniprod, nvmaxprod, jiniprod, jmaxprod
   integer :: iomminprod, iommaxprod, n2prod1, n2prod0
   integer :: nangproj0, nangproj1
   real*8  :: rbalinprod
   real*8  :: Emincut_prod, Emincut_prod_eV

   ! *     for dimension of matrices, total quantites and per processor
   integer :: nomgproc,nanguproc,nanguprocdim,ncanmax
   integer :: nomgprocdim,nprocdim,nomgdim,ncanprocdim
   integer :: ncouprocmax

   real*8 :: rmis1,rfin1,rmis2,rfin2,ah1,ah2,steptot
   real*8, allocatable :: wreal(:), weight(:), cgamma(:), cgamprod(:), weiprod(:)

   integer, allocatable :: ncanproc(:), ibasproc(:, :),&
                          indiomreal(:,:), indangreal(:,:), indproc(:,:),&
                          ntotproc(:),  ncouproc(:), ipcou(:, :)
   integer, allocatable :: iombas(:), nelebas(:)

   ! *     to calculate total memory
   integer*8 :: nointeger_mem,noreal_mem
   integer*8 :: nointegerproc_mem,norealproc_mem

   ! *     number of vibrational states per electronic state

   integer,allocatable :: max_viblevels(:)
   ! ***************************************
   ! *   functions of  mod_gridYpara_01y2  *
   ! ***************************************
   contains
   !--------------------------------------------------
   subroutine input_grid
   !--------------------------------------------------
      implicit none
      integer :: ierror,ir1,ir2,iang,icount,ie
      real*8 :: div
   ! *********************************************************
      namelist /inputgridbase/npun1,rmis1,rfin1,npun1min,                  &
                              npun2,rmis2,rfin2,                           &
                              nangu,                                       &
                              Jtot,iparity,inc,nelecmax,iommin,iommax,j0,  &
                              jini,jmax,nvini,nvmax,                       &
                              nvref,jref,iomref,ielecref
   ! *********************************************************
      namelist /inputprod/iprod,nviniprod,nvmaxprod,max_viblevels,         &
                          jiniprod,jmaxprod,                               &
                          iomminprod,iommaxprod,                           &
                          Rbalinprod,n2prod0,n2prod1,nangproj0,nangproj1,  &
                          Emincut_prod_eV
   ! *********************************************************
      namelist /inputprocess/iphoto,nelec_bnd
   ! *********************************************************
      namelist /inputwrite/iwrt_pot,iwrt_wvp,iwrt_reac_distri,n1plot,n2plot,nangplot
      ! grid step on r1, r2, ang to print

      nointegerproc_mem=0
      norealproc_mem=0
      nointeger_mem=0
      noreal_mem=0

      npun1min=32
      iomref=0
      ielecref=1
      iommin=0

      print *,'(40("="),/,10x,"GridYpara_mod",/,40("="))'
      print *
      print *,'  grid and basis data'
      print *,'  -------------------'
      open(10,file='input.dat',status='old')
      read(10,nml = inputgridbase)
      write(6,nml = inputgridbase)
      flush(6)
      close(10)

      nangu2      = nangu*inc
      iomminprod  = 0
      iommaxprod  = 0
      n2prod1     = npun2
      n2prod0     = 1
      nangproj0   = 1
      nangproj1   = nangu
      iprod       = 0

      allocate(max_viblevels(nelecmax))
      max_viblevels(:)=-1
      print *
      print *,'  products data'
      print *,'  -------------'
      Emincut_prod_eV = 0.d0
      open(10,file='input.dat',status='old')
      read(10,nml = inputprod)

      if(iommaxprod == 0) then
         iommaxprod = min(Jtot, jmaxprod)
      end if

      Emincut_prod = 0.0 !() Emincut_prod_eV/(conve1 * eV2cm) can only be 0.0 as we just assigned Emincut_prod_eV=0.0
      if(nelecmax == 1) then ! QUESTION: What if nelecmax /= 1? 
         max_viblevels(1) = nvmaxprod - nviniprod + 1
      end if

      write(6,nml = inputprod)
      flush(6)
      close(10)

      icount = 0
      do ie = 1, nelecmax
         icount = icount + max_viblevels(ie)
      enddo

      if(icount /= (nvmaxprod - nviniprod + 1)  &
         .and. iprod == 2                       &
         .and. icount > 0) then
         print*, ' in input_grid, namelist /inputprod/ '
         print*, ' (nvmaxprod-nviniprod+1)=', nvmaxprod - nviniprod + 1
         print*, ' no equal to the sum_ie of max_viblevels(ie)'
         
         if(max_viblevels(1) == -1) then ! This only happens if nelecmax /= 1
            print *,' Calculating vibrational states: ', (nvmaxprod-nviniprod+1), ' distributted among all electronic states'
            ! QUESTION: Nothing gets calculated? 
         else
            print *,' correct /inputprod/ and restart '
            stop ' correct /inputprod/ and restart '
         endif
      endif

      nelec_bnd=1
      print *
      print *,'  process data'
      print *,'  ------------'
      open(10,file='input.dat',status='old')
      read(10,nml = inputprocess)
      write(6,nml = inputprocess)
      flush(6)
      close(10)

      print *
      print *,'  write data'
      print *,'  ------------'
      open(10,file='input.dat',status='old')
      read(10,nml = inputwrite)
      write(6,nml = inputwrite)
      flush(6)
      close(10)
      
      if(iwrt_reac_distri == 1)then
         print *, ' iwrt_reac_distri=1 --> writing distriREAC files'
      elseif(iwrt_reac_distri == 2)then
         print *, ' iwrt_reac_distri=2 --> writing Cvj files'
      endif

   ! radial grid integration steps (in angstroms)
      
      if(npun1 > 1)then
         div = dble(npun1 - 1)
         ah1 = (rfin1 - rmis1) / div
      else
         ah1=0.d0
      endif

      if(npun2 > 1) then
         div = dble(npun2 - 1)
         ah2 = (rfin2 - rmis2) / div
      else
         ah2 = 0.d0
      endif

      if(npun1 > 1 .and. npun2 > 1) then
         steptot = ah1 * ah2
      elseif(npun1 > 1  .and. npun2 <= 1)then
         steptot=ah1
      elseif(npun1 <= 1 .and. npun2 > 1)then
         steptot=ah2
      else
         steptot=1.d0
      endif

!     angular grid: Gauss legendre

      allocate(wreal(nangu2),weight(nangu2),cgamma(nangu2), &
               cgamprod(nangu2),weiprod(nangu2),            &
               stat=ierror)

      norealproc_mem = norealproc_mem + 5 * nangu2
      print *, 'norealproc_mem= ',norealproc_mem

      if(nangu == 1) then
         wreal(1)=1.d0
         weight(1)=1.d0
         cgamma(1)=1.d0
      else 
         if (idproc == 0) then
            print *, "inc = ", inc, "j0 = ", j0
         end if

         call gauleg(weight,cgamma,nangu2)

         do iang=1, nangu2
            cgamprod(iang) = cgamma(iang)
            weiprod(iang) = weight(iang)
            wreal(iang) = weight(iang) * inc
         enddo

         if(idproc == 0) then
            print("(//,10x,'Angular grid to plot pes and wvp every ',i3,//)"),nangplot
            do iang=1, nangu, nangplot
               print *, iang,dacos(cgamma(iang))*180.d0/pi
            enddo
         endif
      endif

!     basis set conditions

      if(abs(iparity) /= 1) then
         print *,'ipar= ',iparity
         print *,'|iparity| must be 1 '
         stop
      endif
      
      ! la primera condicion de esta comparación no parece tener sentido (nangu2 = nangu * inc) -> (nangu*inc-1)*inc ??
      if(jmax > (nangu2-1)*inc .and. npun2 > 1) then
         print *, "It must be true that jmax < (nangu2-1)*inc, but:"
         print *, 'jmax= ',jmax,' > (nangu2-1)*inc= ',(nangu2-1)*inc
         stop
      endif

      return
   end subroutine input_grid

!--------------------------------------------------
!--------------------------------------------------

end module mod_gridYpara_01y2
