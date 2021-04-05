subroutine conservation_compute_energies
! this routine computes the various energies to check
! for energy conservation
  
  use blmod, only: nt, delta_mass, cmass, cr, gravity_switch, eps, vel, r, p, &
                    total_initial_energy, time, tdump_scalar, lum_observed, &
                    dtime, energy_from_heating, radiated_energy, simple_heating,&
                    pdVwork
  use parameters
  use physical_constants
  implicit none

  integer :: i
  logical :: force
  logical :: outputflag

  real*8 :: egrav, eint, ekin

!------------------------------------------------------------------------------
  
  egrav = 0.0d0
  eint = 0.0d0
  ekin = 0.0d0

  do i=1,imax
     egrav = egrav - ggrav*delta_mass(i)*(cmass(i)+mass_gravity_switch*mass_extragravity*msun)/cr(i) *gravity_switch
     eint = eint + eps(i)*delta_mass(i)
  enddo


  
  do i=1,imax-1
     ekin = ekin + 0.5d0*(0.50d0*(vel(i+1)+vel(i)))**2 * delta_mass(i)
  enddo
  ekin = ekin + 0.5d0*(vel(imax))**2 * delta_mass(imax)

  radiated_energy = radiated_energy + dtime*lum_observed
  energy_from_heating = energy_from_heating + dtime*sum(simple_heating(1:imax-1)*delta_mass(1:imax-1))
  !pdVwork = pdVwork + p(imax)*4*pi*r(imax)**2*vel(imax)*dtime
  pdVwork = pdVwork + p(1)*4*pi*r(1)**2*vel(1)*dtime

  if(time.eq.0.0d0) then
      total_initial_energy = egrav+eint+ekin
  endif

!#if 0
!  if(mod(nt,1000).eq.0.or.force) then
!     write(6,"(A18,A18,A18,A18)") "egrav","eint","ekin","etot"
!     write(6,"(1P10E18.9)") egrav,eint,ekin,egrav+eint+ekin
!  endif
!#endif


    outputflag = .false.

    if(mod(nt,10000) .eq. 0) then
       outputflag = .true.
    endif

    if(time.eq.0.0d0.or.time.gt.tdump_scalar) then
       outputflag = .true.
    endif



    if(time .eq. 0.0d0) then
      open(666,file=trim(adjustl(outdir))//"/conservation.dat",&
              status='unknown',position='append')
      write(666,*) 'time,egrav,eint,ekin,egrav+eint+ekin,pdVwork, &
              egrav+eint+ekin-total_initial_energy, radiated_energy, energy_from_heating'
      close(666)

    end if

    if (outputflag) then
        open(666,file=trim(adjustl(outdir))//"/conservation.dat",&
                status='unknown',position='append')
        write(666,"(1P10E18.9)") time,egrav,eint,ekin,egrav+eint+ekin,pdVwork, &
                egrav+eint+ekin-total_initial_energy, radiated_energy, energy_from_heating
        close(666)
    end if









end subroutine conservation_compute_energies

