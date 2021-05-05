subroutine conservation_compute_energies
! this routine computes the various energies to check
! for energy conservation
  
  use blmod, only: nt, delta_mass, cmass, cr, gravity_switch, eps, vel, r, p, &
                    total_initial_energy, time, tdump_scalar, lum_observed, &
                    dtime, energy_from_heating, radiated_energy, simple_heating,&
                    pdVwork_inner,pdVwork_outer
  use parameters
  use physical_constants
  implicit none

  integer :: i
  logical :: force
  logical :: outputflag

  real*8 :: egrav, eint, ekin, pdVterm, E1, E2

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
  pdVwork_outer = pdVwork_outer + p(imax)*4*pi*r(imax)**2*vel(imax)*dtime
  pdVwork_inner = pdVwork_inner + p(1)*4*pi*r(1)**2*vel(1)*dtime
  pdVterm = pdVwork_inner - pdVwork_outer !energy put into the ejecta due to pressure work at boundaries
  
  E1 = egrav+eint+ekin  !total energy of the ejecta
  
  if(time.eq.0.0d0) then
      total_initial_energy = E1
  endif

  E2 = total_initial_energy + pdVterm + energy_from_heating - radiated_energy
  !the initial ejecta energy plus energy input and output up to t=time,
  !theoretically E2 = E1 if energy is conserved.

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
      write(666,*) 'time,egrav,eint,ekin,E1,pdVwork_outer,pdVwork_inner, &
              pdVterm,energy_from_heating,radiated_energy,E2,E1-E2'
      write(666,*) 'E1 = egrav+eint+ekin is the total ejecta energy, E2 = initial total ejecta energy &
              + pdVterm + energy_from_heating - radiated_energy, E1=E2 if energy is conserved.'
      close(666)

    end if

    if (outputflag) then
        open(666,file=trim(adjustl(outdir))//"/conservation.dat",&
                status='unknown',position='append')
        write(666,"(1P15E18.9)") time,egrav,eint,ekin,E1,pdVwork_outer, &
                pdVwork_inner, pdVterm, energy_from_heating, radiated_energy, E2, E1-E2
        close(666)
    end if









end subroutine conservation_compute_energies

