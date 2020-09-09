## Homologous expansion

###Changes

```fortran
parameters:
+ bomb_mode=2
bomb_tend           = 0.0d0
bomb_mass_spread    = 0.0d0 #(in solar mass)

final_energy=0.0d0
boxcar_smoothing=0
profile_name 		= "profiles/Homologous_expansion.dat"
comp_profile_name	= "profiles/Homologous_expansion_composition.dat"

eoskey = 1
radiation = 0
Ni_switch = 0
 + Ni_by_hand = 0
mass_excision = 0



```



```fortran
problem.F90
!  call compose_opacity_tables_OPAL
!
!  call opacity(rho(:),temp(:),kappa(:),kappa_table(:),dkappadt(:))
!
!  call optical_depth(rho(:), r(:), kappa_table(:), tau(:))
!
!  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))
!
!
!  call read_BolCorr
```

```fortran
hydro.F90
 !calculate heating term due to Ni
!  if(time.ge.time_Ni) then
!      time_Ni = time_Ni + Ni_period
!      call nickel_heating
!  endif
Ni_heating(:)=0.0d0
  !calculate heating term due to bomb
!  if(do_bomb .and. time.ge.bomb_tstart .and. time.le.bomb_tend) then
!      call bomb_pattern
!  else
!      bomb_heating(:) = 0.0d0
!  endif




!  call opacity(rho(:),temp_temp(:),kappa(:),kappa_table(:),dkappadt(:))
!
!  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))
```

```fortran
blstep.F90

!     if(.not.sedov) then
!        call analysis
!     endif
```





```fortran
output.F90
! filename = trim(adjustl(outdir))//"/Ni_deposit_function.xg"
! call output_single_mass(Ni_deposit_function,filename)

! filename = trim(adjustl(outdir))//"/He_1.xg"
! call output_single_mass(ion_fractions(He_number,1,:),filename)

! filename = trim(adjustl(outdir))//"/He_2.xg"
! call output_single_mass(ion_fractions(He_number,2,:),filename)

! filename = trim(adjustl(outdir))//"/He_3.xg"
! call output_single_mass(ion_fractions(He_number,3,:),filename)

! filename = trim(adjustl(outdir))//"/H_1.xg"
! call output_single_mass(ion_fractions(H_number,1,:),filename)

! filename = trim(adjustl(outdir))//"/H_2.xg"
! call output_single_mass(ion_fractions(H_number,2,:),filename)

! filename = trim(adjustl(outdir))//"/free_electron_frac.xg"
! call output_single_mass(free_electron_frac,filename)



! filename = trim(adjustl(outdir))//"/photosphere_tracer.xg"
    ! call output_single_mass(photosphere_tracer,filename)
```