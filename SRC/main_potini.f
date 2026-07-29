      program potini
!     *********************************************************************
!     Initialize the potential, reactant and product wf, and dipole transition moment
!     to be read to proceed for several wave packet calculations witj MadWave3
!     for AB+C --> AC + C
!     for 01+2 --> 02 + 1
!     for different initial states
!     *********************************************************************

      use mod_gridYpara_01y2
      use mod_pot_01y2
      use mod_baseYfunciones_01y2
      use mod_photoini_01y2
      use mod_Hphi_01y2
      use mpi
      
      implicit none
      character*40 filename
      integer ierr,ielec


      
      ! Initialize MPI environment and get proc's ID and number of proc in
      ! the partition.

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

      write(filename,'("salpot.",i3.3)')idproc
      open(6,file=filename,status='unknown')
      print '(40("="),//)'
      print '(10x,"Potini for Madwave3 version 6 ",//)'
      print *, ' output of proc. idproc= ',idproc,' of nproc= ',nproc
      print '(/,40("="),//)'
      
!     initialization of data

      call input_grid()
      call pot0()
      call paralelizacion()
      print *, 
      print *, '  --- Initialization of the potential ---'
      print *, 

      call setxbcpotele(iomdiat,iomatom,sigdiat,sigatom,nelec,nelecmax)
      print *, 
      print *, '  --- end initialization of the potential ---'
      print *, 
      if(nelec.gt.nelecmax.or.nelec.eq.0)then
          print *,'  !!! nelec= ',nelec
     &                            ,'  while nelecmax= ',nelecmax
          call flush(6)
          stop
      endif
         open(10,file='pot/cont.pot',status='new')
         write(10,*)nelec
         do ielec=1,nelec
            write(10,*)iomdiat(ielec),iomatom(ielec)
     &           ,sigdiat(ielec),sigatom(ielec)
         enddo
         close(10)

!     determining basis

      call basis
      
! reactants functions calculation

      call angular_functions

      if(npun1.eq.1)then

         write(6,*)' --> calling rigid rotor energies <-- '
         call flush(6)
         call radial_functions01eq_write

      else
         call radial_functions01_write

         if(iprod.eq.2)then
 
            call product_radialf_write

         endif
      endif

!     determining potential

      call pot1

!     determining electric dipole transition

      if(iphoto.gt.0.and.iphoto.le.2)then
         call pot2                ! to initialize global indexes
         if(iphoto.eq.1)then
            call set_trans_dipole
         endif
      elseif(iphoto.eq.3)then
         call pot2                ! to initialize global indexes
         call set_coupling        ! for electronic predissociation
      endif
      
      stop
      end program
