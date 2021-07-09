subroutine get_inner_outer_mass_from_profile(prof_name,inner_mass,outer_mass)

  implicit none

!Input:
  character(*)         :: prof_name

!Output:
  real*8               :: inner_mass
  real*8               :: outer_mass

!Local:
  integer              :: profile_zones
  integer              :: i,ibuffer
  real*8               :: buffer
  real*8, allocatable  :: pmass(:)

!------------------------------------------------------------------------------

  open(666,file=trim(prof_name),status='unknown',&
       form='formatted',action='read')
  read(666,*)
  read(666,*) profile_zones
  allocate(pmass(profile_zones))
  
  do i=1,profile_zones
    read(666,*) ibuffer,pmass(i)
  enddo
  close(666)

  outer_mass = pmass(profile_zones)
  inner_mass = pmass(1)

  deallocate(pmass)

end subroutine get_inner_outer_mass_from_profile

!******************************************************************************

subroutine get_inner_outer_radius_from_profile(prof_name,inner_mass, &
        inner_radius,outer_radius)

  implicit none

!Input:
  character(*)         :: prof_name
  real*8               :: inner_mass

!Output:
  real*8               :: inner_radius
  real*8               :: outer_radius

!Local:
  integer              :: profile_zones
  integer              :: i,ibuffer
  real*8               :: buffer
  real*8, allocatable  :: pradius(:)
  real*8, allocatable  :: pmass(:)

!------------------------------------------------------------------------------

  open(666,file=trim(prof_name),status='unknown',&
       form='formatted',action='read')
  read(666,*)
  read(666,*) profile_zones
  allocate(pradius(profile_zones))
  allocate(pmass(profile_zones))

  do i=1,profile_zones
     read(666,*) ibuffer, pmass(i), pradius(i)
  enddo
  close(666)
  
  !inner mass takes into account the excised mass, and does not necessarily
  !coincide with the inner mass of the profile
  call map_map(inner_radius,inner_mass,pradius,pmass,profile_zones)

  outer_radius = pradius(profile_zones)
  
  deallocate(pradius)
  deallocate(pmass)

end subroutine get_inner_outer_radius_from_profile

!******************************************************************************

subroutine get_ncomps_from_profile(prof_name,xncomps)

  implicit none

!Input:
  character(*)         :: prof_name

!Output:
  integer              :: xncomps

!Local:
  integer              :: profile_zones

!------------------------------------------------------------------------------

  open(666,file=trim(prof_name),status='unknown',&
       form='formatted',action='read')
  read(666,*) profile_zones,xncomps

  close(666)

end subroutine get_ncomps_from_profile

!****************************** READ PROFILE **********************************
      
subroutine read_profile(prof_name)

  use blmod, only: mass, cmass, vel, rho, temp, ncomps, ye,ye_initial,abar, comp_details,&
                    eps, p, cs2, dedt, dpdt, entropy, zav, p_rad, entropy_frominput,&
                    expansion_timescale, vel_inner
  use parameters
  use physical_constants
  use eosmodule, only: init_ionpot
  implicit none

  character(*) :: prof_name
  integer :: profile_zones
  integer :: i,l
  integer :: ibuffer
  integer :: keytemp, keyerr

  real*8,allocatable :: pmass(:), pradius(:), ptemp(:), prho(:), pvel(:),pye(:),pentropy(:),pexpansion_timescale(:)

!------------------------------------------------------------------------------

  open(666,file=trim(prof_name),status='unknown',form='formatted',action='read')
  read(666,*)
  read(666,*) profile_zones
  profile_zones = profile_zones
  write(*,*) "We have ",profile_zones, "profile zones."

  !----------------- read the profile and map it on the grid --------------------
  allocate(pmass(profile_zones))
  allocate(pradius(profile_zones))
  allocate(ptemp(profile_zones))
  allocate(prho(profile_zones))
  allocate(pvel(profile_zones))
  allocate(pye(profile_zones))
  allocate(pentropy(profile_zones))
  allocate(pexpansion_timescale(profile_zones))

  do i=1,profile_zones
     read(666,*) ibuffer, pmass(i), pradius(i), ptemp(i), prho(i), pvel(i), pye(i), pentropy(i), pexpansion_timescale(i)
  enddo

  do i=1,imax !velocity lives at the cell edges
      call map_map(vel(i), mass(i),pvel,   pmass,profile_zones)
  enddo
  vel_inner = vel(1)

  do i=1,imax-1 !temperature and density live at the cell centers
     call map_map(rho(i), cmass(i),prho,   pmass,profile_zones)
     call map_map(temp(i),cmass(i),ptemp,  pmass,profile_zones)
     call map_map(entropy_frominput(i),cmass(i),pentropy,  pmass,profile_zones)
     call map_map(expansion_timescale(i),cmass(i),pexpansion_timescale, pmass,profile_zones)
  enddo

  entropy_frominput(imax) = entropy_frominput(imax-1)
  expansion_timescale(imax) = expansion_timescale(imax-1)

  do i=1,imax
    expansion_timescale(i) = expansion_timescale(i)*1000.0d0 ! s to ms
  enddo

  if(continuous_boundary_switch)then
    rho(imax)=rho(imax-1)!passive boundary condition
    temp(imax)=temp(imax-1)!passive boundary condition
  else
    rho(imax) = 0.0d0 !passive boundary condition
    temp(imax) = 0.0d0 !passive boundary condition
  end if


!!!------------------------- read composition profile ---------------------------
  if(ncomps.gt.0) then
     call read_profile_compositions(composition_profile_name)

     if(eoskey.eq.2) then
       ! initialize some variables need in the
       ! saha solver -- need to have composition info at this point
       call init_ionpot
     endif

     !initialize zav (for the Saha solver, Paczynski EOS), assuming full ionization
     do l=1, ncomps
       do i=1, imax
         zav(l,i) = comp_details(l,2)
       end do
     end do
  endif
!  for kilonovae, no composition profile
  if(ncomps.eq.0) then
    do i = 1,imax-1
      !ye lives at the cell centers
      call map_map(ye_initial(i), cmass(i),pye,pmass,profile_zones)
    end do
    ye_initial(imax) = ye_initial(imax-1)

    abar(:) = mu  !print*,"Assume mean molecular weight = 150!"
    ye(:) = ye_afternucleosynthesis

  end if

  deallocate(pmass)
  deallocate(pradius)
  deallocate(ptemp)
  deallocate(prho)
  deallocate(pvel)
  deallocate(pye)
  deallocate(pentropy)
  deallocate(pexpansion_timescale)


!------------- find other hydrodynamical quantities from the EOS --------------


  !call equation of state
  keytemp = 1

  call eos(rho(1:imax-1),temp(1:imax-1),ye(1:imax-1), &
         abar(1:imax-1),p(1:imax-1),eps(1:imax-1), &
         cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), &
         entropy(1:imax-1),p_rad(1:imax-1),keyerr,keytemp,eoskey)


  if(continuous_boundary_switch) then
    eps(imax)=eps(imax-1)!passive boundary condition
    p(imax)=p(imax-1)!active boundary condition, used in the velocity update
  else
    eps(imax) = 0.0d0!passive boundary condition
    p(imax) = 0.0d0!active boundary condition, used in the velocity update
  end if


end subroutine read_profile


!******************************************************************************
subroutine map_linterp(x1,x2,y1,y2,x,y)

! perform linear interpolation      
  implicit none

  real*8 :: slope,x1,x2,y1,y2,x,y

    if (x2.lt.x1) then
       stop "Error in linterp!"
    endif

    if (x2.ne.x1) then
        slope = (y2 - y1) / (x2 - x1)
    else
        slope = 0
    endif

    y = slope*(x-x1) + y1

end subroutine  map_linterp

!******************************************************************************
subroutine map_find_index(zones,array,goal,upper_index,lower_index)

  ! bisection search
  implicit none

  integer :: zones,i
  real*8 :: array(*)
  real*8 :: goal
  integer :: middle_index,upper_index,lower_index

    if(goal.lt.array(1) .or. goal.gt.array(zones)) then
        write(*,*) 'value passed to map_find_index is out of the array'
        stop
    end if

    lower_index = 1
    upper_index = zones

    do while ( (upper_index - lower_index) .gt. 1 )
       middle_index = (lower_index + upper_index) * 0.5d0
       if ( goal .le. array(middle_index) ) then
            upper_index = middle_index
       else
            lower_index = middle_index
       end if
    enddo

end subroutine map_find_index


!******************************************************************************

subroutine map_map(point_value,point_radius,parray,pradius,zones)

  implicit none

  real*8 :: point_value, point_radius
  real*8 :: pradius(*), parray(*)
  integer :: zones
  integer :: upper_index, lower_index

    if (point_radius .ge. pradius(1) .and. &
       point_radius .le. pradius(zones) )  then

     call map_find_index(zones,pradius,point_radius, &
          upper_index,lower_index)

     call map_linterp( pradius(lower_index),pradius(upper_index), &
          parray(lower_index), parray(upper_index),  & 
          point_radius, point_value )

    else if (point_radius .lt. pradius(1)) then
     ! linear extrapolation
     call map_linterp(pradius(1),pradius(2), & 
          parray(1),parray(2),point_radius,point_value)

    else if (point_radius .gt. pradius(zones)) then
     ! linear extrapolation
     call map_linterp(pradius(zones-1),pradius(zones), & 
          parray(zones-1),parray(zones),point_radius,point_value)
    endif


end subroutine map_map

!******************************************************************************
subroutine integrate_radius_initial

  use blmod, only: r, rho, cr, delta_mass
  use parameters
  use physical_constants
  implicit none

  integer :: i

!------------------------------------------------------------------------------
  
  do i=1,imax-1
   r(i+1) = ( 3.0d0/(4.0d0*pi) * delta_mass(i)/rho(i) + r(i)**3 )**(1.0d0/3.0d0)



   if(rho(i).le.0.0d0) then
      stop "Negative density. Profile extends not to large enough radii."
   endif
  enddo


  do i=1,imax-1
    cr(i) = ( ( r(i)**3 + r(i+1)**3 ) / 2.0d0 )**(1.0d0/3.0d0)
  enddo

  cr(imax) = r(imax) + (r(imax) - cr(imax-1))
  !passive boundary condition, used in the expression for the velocity update,
  !but multiplied by the artificial viscosity, which is zero at the last point



end subroutine integrate_radius_initial
