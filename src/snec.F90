program snec

  use blmod, only: dtime, dtime_p, time, nt, ntstart, tstart,   &
     tdump, tdump_scalar, rho, tdump_check,lum_photo,eps,ye,mass,index_photo
  use parameters
  use outinfomod, only: outinfo_count
  use heating_rate_LR15_module

  implicit none
  external :: Blackbody
  logical :: OutputFlag = .false.
  logical :: OutputFlagScalar = .false.
  logical :: OutputFlagCheck = .false.
!  integer :: test_maxcount = 5
!  integer :: test_count = 0

  integer :: system_time(8)
  real*8 :: extra_output_points(51)
  integer :: i



  !------------------------------------------------------------------------------

  call date_and_time(values=system_time)
  print '(i4, 5(a, i2.2), " UTC")', system_time(1), '/', system_time(2), '/', system_time(3), ' ', &
          system_time(5), ':', system_time(6), ':', system_time(7)


  write(*,*) "***********************************"
  write(*,*) "* Supernova Explosion Code (SNEC) *"
  write(*,*) "***********************************"

  write(*,*)

  do i =1,51
      extra_output_points(i) = 24*3600*10**(-7.0d0+(i-1)*0.1d0)
  end do

  i = 1

  ! *****************************************************
! INITIALIZATION
! *****************************************************        


  call input_parser

  call problem

  call artificial_viscosity
  
! output before first timestep
  call output_all(0)
  call output_all(1)
  call output_all(2)
  
  call timestep

  tdump_check = tstart+dtout_check
  tdump_scalar = tstart+dtout_scalar
  tdump = tstart+dtout
  time = tstart
  nt = ntstart

!  ye(:) = 62.0/150.0
!  print*,"ye for EOS is set to be 62/150!"
! *****************************************************
! MAIN LOOP
! *****************************************************

  IntegrationLoop: do
!      print*,'time:',time
!     test_count = test_count +1
!     if(test_count>test_maxcount) then
!         exit
!     end if

     dtime_p = dtime
     ! determine dt
     call timestep


     if(ntinfo.gt.0) then
        if(mod(nt,ntinfo).eq.0) then
           ! print useful info to stdout
           call outinfo
        endif
     endif

     if((time+dtime).gt.tend) dtime = tend-time

     ! actual integration step
     call blstep

     ! increment timestep
     nt = nt + 1

     ! various output related things
     if (ntout.gt.0) then
        if ( mod(nt,ntout) .eq. 0) OutputFlag = .true.
     endif

     if (ntout_scalar.gt.0) then
        if ( mod(nt,ntout_scalar) .eq. 0 ) OutputFlagScalar = .true.
     endif

     if (ntout_check.gt.0) then
        if ( mod(nt,ntout_check) .eq. 0 ) OutputFlagCheck = .true.
     endif

     if ( time.ge.tdump) then
        tdump=tdump+dtout
        OutputFlag = .true.
     endif

     if ( time.ge.tdump_scalar) then
        tdump_scalar=tdump_scalar+dtout_scalar
        OutputFlagScalar = .true.
     endif

     if ( time.ge.tdump_check) then
        tdump_check=tdump_check+dtout_check
        OutputFlagCheck = .true.
     endif

      !The following is used to add output at very early times, extra_output_points are
      ! spaced between 10^(-7) day to 10^(-2) day on a log scale.
     if (i .le. 51) then
         if (time .ge. extra_output_points(i))  then
             OutputFlagScalar = .true.
             i = i + 1
         end if
     end if
     ! Besides, we need to ensure the location of photosphere is fully resolved
     !so that the light curves can be trusted.
     !So we add another threshold for output: index_photo < imax - 5
     if(index_photo .ge. imax-5) then
         OutputFlagScalar = .false.
     end if



     ! increment time
     time = time+dtime

     if (OutputFlag) then
        call output_all(0)
        call output_all(1)
        OutputFlag = .false.
     endif

     if (OutputFlagScalar) then
        call output_all(2)
        OutputFlagScalar = .false.
     endif

     if (OutputFlagCheck) then
        OutputFlagCheck = .false.
     endif

     if (time.eq.tend) then
        write(*,*) "Done! :-) tend reached"
        call output_all(0)
        call output_all(1)
        call output_all(2)
        exit
     endif

     if (nt.ge.ntmax) then
        write(*,*) "Done! :-) ntmax reached"
        call output_all(0)
        call output_all(1)
        call output_all(2)
        exit
     endif



  enddo IntegrationLoop


  call date_and_time(values=system_time)
  print '(i4, 5(a, i2.2), " UTC")', system_time(1), '/', system_time(2), '/', system_time(3), ' ', &
          system_time(5), ':', system_time(6), ':', system_time(7)


end program snec
