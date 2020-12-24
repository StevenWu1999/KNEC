subroutine analysis

  use blmod
  use parameters
  use physical_constants
  implicit none

  character(len=1024) :: filename

  integer :: i,j
  integer :: lum_observed_maxindex
  real*8 :: ratio,flux
  real*8 :: ratio_color(5)
  real*8 :: ratio_color_mag(5)
  real*8 :: flux_color(5)


  real*8 :: AB_zeropoint(5)
!  real*8 :: T_eff_for_BC
!  real*8 :: bol_corr_used(11)
!  real*8, parameter :: T_eff_min = 5000.0d0
  real*8 :: lum_i
  real*8 :: D !Distance between source and observer, 40 Mpc by default

  D = 40 *Mpc
!------------------------------------------------------------------------------
!---------- Calculating optical depth and tracing the photosphere -------------


  !kappa_table (without the opacity floor) is used to trace the photosphere
  call optical_depth(rho(:), r(:), kappa_table(:), tau(:))

  !find the grid point, where the photosphere is located
  do i=imax-1, 1, -1
     if(tau(i).gt.0.66d0) then
        index_photo = i + 1
        exit
     endif
  enddo

  !fix the moment, when the photosphere reaches the inner boundary, if it does
  if(tau(1).lt.0.66d0) then
     index_photo = 1
     if(photosphere_fell_on_the_center.eq.0) then
        photosphere_fell_on_the_center = 1
        open(unit=666, &
             file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
             status="unknown",form='formatted',position="append")
        write(666,*) 'Photosphere reached the center at ', time, 'seconds'
        close(666)
     end if
  end if

  !find the values of some variables at the photosphere
  if(photosphere_fell_on_the_center.eq.0) then
     call map_map(lum_photo,  0.66d0, lum(imax:1:-1),  tau(imax:1:-1),imax)
     call map_map(mass_photo, 0.66d0, mass(imax:1:-1), tau(imax:1:-1),imax)
     call map_map(vel_photo,  0.66d0, vel(imax:1:-1),  tau(imax:1:-1),imax)
     call map_map(rad_photo,  0.66d0, r(imax:1:-1),    tau(imax:1:-1),imax)
  else
     lum_photo   = lum(1)
     mass_photo  = mass(1)
     vel_photo   = vel(1)
     rad_photo   = r(1)
  end if


  rad_max = r(imax)
  !photosphere tracer is equal to 1 at the photosphere, and 0 everywhere else
  !used for visualization
  photosphere_tracer(:) = 0.0d0
  photosphere_tracer(index_photo) = 1.0d0

  !check if the photosphere moves through the regions with wrong opacity
  !log10(T) = 3.75 is the lower boundary of the OPAL opacity tables
!  if(metallicity(index_photo).gt.envelope_metallicity .and. &
!       log10(temp(index_photo)).lt.3.75d0 .and. &
!       kappa_table(index_photo).gt.kappa(index_photo) &
!       .and. index_photo.gt.1 ) then
!     opacity_corrupted = 1
!  else
!     opacity_corrupted = 0
!  end if


!------------------------ Tracing the luminosity shell ------------------------
  index_lumshell = imax
  do i=1, imax - 1
     if(tau(imax-i).gt.(clite/vel(imax-i))) then
        index_lumshell = imax - i + 1
        exit
     end if
  end do
  if (tau(2).le.(clite/vel(2))) then
     index_lumshell = 1
  end if
  mass_lumshell = mass(index_lumshell)

  !characteristic diffusion and expansion times for different shells
  do i=1, imax-1
     time_diff(i) = kappa(i)*rho(i)*(r(imax)-r(i))**2.0/clite
     time_exp(i) = (r(imax)-r(i))/(vel(imax)-vel(i))
  end do

  !internal energies of shells from the given radius out to the surface
  do i=1, imax
     E_shell(i) = sum(eps(i:imax)*delta_mass(i:imax))
  end do



!------------------- Calculate the observed luminosity ------------------------

!  observed luminosity is the sum of lum_photosphere and Ni contribution
!  if(photosphere_fell_on_the_center.eq.0) then
!     lum_observed = lum_photo + sum(Ni_energy_rate* &
!          Ni_deposit_function(index_photo:imax)*delta_mass(index_photo:imax))
!  else
!     lum_observed = lum(1) + &
!          sum(Ni_energy_rate*Ni_deposit_function(1:imax)*delta_mass(1:imax))
!  end if


  lum_observed_maxindex = 1000
  if(photosphere_fell_on_the_center.eq.0) then
      lum_observed = lum_photo + sum(simple_heating(index_photo:imax)*delta_mass(index_photo:imax))
      if(index_photo<lum_observed_maxindex) then
          lum_observed_min = lum_photo + &
                  sum(simple_heating(index_photo:lum_observed_maxindex)*delta_mass(index_photo:lum_observed_maxindex))
      else
          lum_observed_min = lum_photo + sum(simple_heating(index_photo:imax)*delta_mass(index_photo:imax))
      end if

  else
      lum_observed = lum(1) + &
              sum(simple_heating(1:imax)*delta_mass(1:imax))
      if(index_photo<lum_observed_maxindex) then
          lum_observed_min = lum(1) + &
                  sum(simple_heating(index_photo:lum_observed_maxindex)*delta_mass(index_photo:lum_observed_maxindex))
      else
          lum_observed_min = lum(1) + sum(simple_heating(index_photo:imax)*delta_mass(index_photo:imax))
      end if

  end if




  !write down the time when the contribution of the Ni above the
  !photosphere to the luminosity is greater than 5%
!  if(shockpos_stop.eq.1 .and. &
!       abs((lum_observed - lum_photo)/lum_photo).gt.0.05 .and. &
!       Ni_contributes_five_percents.eq.0) then
!
!     Ni_contributes_five_percents = 1
!     open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
!          status="unknown",form='formatted',position="append")
!     write(666,*) 'Ni contribution to the luminosity is 5% at ', time, 'seconds'
!     close(666)
!
!  end if

!-------- Calculation of color magnitudes using bolometric corrections --------

  T_eff = (lum_photo/(4.0d0*pi*sigma_SB*rad_photo**2))**0.25d0

  !see Eq.(3) of Swartz et al., ApJ 374:266 (1991) and explanation there
!  T_eff_for_BC = MAX(T_eff,T_eff_min)

  do i=1,index_photo-1
      Temp_for_color(i) = T_eff
  end do

  do i=index_photo,imax-1
      lum_i = delta_mass(i)*simple_heating(i)
      Temp_for_color(i) = (lum_i/(4.0d0*pi*sigma_SB*r(i)**2))**0.25d0
  end do


  !here in cases, when the effective temperature goes beyond the boundaries
  !of the table BolCorr.dat, the linear extrapolation is used
!  do i=1, 11
!     call map_map(bol_corr_used(i),T_eff_for_BC,bol_corr(1:nlines_bol_corr,i), &
!          temp_bol_corr,nlines_bol_corr)
!     magnitudes(i) = sun_mag - bol_corr_used(i) - 2.5d0*log10(lum_photo/sun_lum)
!  end do

  !write down the magnitudes in a file
!  if(time.eq.0.0d0.or.time.gt.tdump_scalar) then
!     open(666,file=trim(adjustl(outdir))//"/magnitudes.dat",&
!          status='unknown',position='append')
!     write(666,"(13E18.9)") time, T_eff_for_BC, magnitudes(1:11)
!     close(666)
!  endif




  if(time.eq.0.0d0.or.time.gt.tdump_scalar) then

      ! color light curves
    !   lum_color(:) = 0.0d0
    !   flux_color(:) = 0.0d0
    !   call Filter(T_eff,ratio_color,ratio_color_mag,AB_zeropoint)
    !   lum_color = lum_color + lum_photo*ratio_color
    !   flux_color = flux_color + lum_photo*ratio_color_mag/(4*pi*(40*Mpc)**2)

    !   do i = index_photo,imax-1
    !       call Filter(temp(i),ratio_color,ratio_color_mag,AB_zeropoint)
    !       lum_color = lum_color + simple_heating(i)*delta_mass(i)*ratio_color
    !       flux_color = flux_color + simple_heating(i)*delta_mass(i)*ratio_color_mag/(4*pi*(40*Mpc)**2)
    !   end do

    !   Magnitude_UBVRI = -2.5d0* log10(flux_color/AB_zeropoint)

    !   open(666,file=trim(adjustl(outdir))//"/lum_color.dat",&
    !           status='unknown',position='append')
    !   write(666,"(13E18.9)") time, lum_color(1:5)
    !   close(666)


    !   open(666,file=trim(adjustl(outdir))//"/Magnitude_UBVRI.dat",&
    !           status='unknown',position='append')
    !   write(666,"(13E18.9)") time, Magnitude_UBVRI(1:5)
    !   close(666)

      !Gemini_ugriz
      do j = 1,5
          flux = 0.0d0

          if (photosphere_fell_on_the_center .eq. 0) then
              ! Blackbody assumption for the photosphere:
              call Blackbody_lambda_integral(T_eff,Wavelength_Gemini_ugriz(:,j), Trans_Gemini_ugriz(:,j),&
                      Minmax_Gemini_ugriz(:,j),size(Wavelength_Gemini_ugriz(:,j)),ratio)
              flux = flux + ratio*lum_photo/(4*D**2) ! not 4 pi D^2 because pi has been reduced

          end if

          ! Blackbody assumption for each radioactive layer above the photosphere:
          do i = index_photo,imax-1
              call Blackbody_lambda_integral(temp(i),Wavelength_Gemini_ugriz(:,j), Trans_Gemini_ugriz(:,j),&
                      Minmax_Gemini_ugriz(:,j),size(Wavelength_Gemini_ugriz(:,j)),ratio)
              flux = flux + ratio*simple_heating(i)*delta_mass(i)/(4*D**2)
          end do

          Magnitude_Gemini_ugriz(j) = -2.5d0*log10(flux/ABzeropoint_Gemini_ugriz(j))

      end do


      !Gemini_JHKs
      do j = 1,3
          flux = 0.0d0
          if (photosphere_fell_on_the_center .eq. 0) then
              ! Blackbody assumption for the photosphere:
              call Blackbody_lambda_integral(T_eff,Wavelength_Gemini_JHKs(:,j), Trans_Gemini_JHKs(:,j),&
                      Minmax_Gemini_JHKs(:,j),size(Wavelength_Gemini_JHKs(:,j)),ratio)
              flux = flux + ratio*lum_photo/(4*D**2) ! not 4 pi D^2 because pi has been reduced
          end if

          ! Blackbody assumption for each radioactive layer above the photosphere:
          do i = index_photo,imax-1
              call Blackbody_lambda_integral(temp(i),Wavelength_Gemini_JHKs(:,j), Trans_Gemini_JHKs(:,j),&
                      Minmax_Gemini_JHKs(:,j),size(Wavelength_Gemini_JHKs(:,j)),ratio)
              flux = flux + ratio*simple_heating(i)*delta_mass(i)/(4*D**2)
          end do
          Magnitude_Gemini_JHKs(j) = -2.5d0*log10(flux/ABzeropoint_Gemini_JHKs(j))

      end do


    !CTIO_BVRIJHK
      do j = 1,7
        flux = 0.0d0
        if (photosphere_fell_on_the_center .eq. 0) then
            ! Blackbody assumption for the photosphere:
            call Blackbody_lambda_integral(T_eff,Wavelength_CTIO_BVRIJHK(:,j), Trans_CTIO_BVRIJHK(:,j),&
                    Minmax_CTIO_BVRIJHK(:,j),size(Wavelength_CTIO_BVRIJHK(:,j)),ratio)
            flux = flux + ratio*lum_photo/(4*D**2) ! not 4 pi D^2 because pi has been reduced
        end if

        ! Blackbody assumption for each radioactive layer above the photosphere:
        do i = index_photo,imax-1
            call Blackbody_lambda_integral(temp(i),Wavelength_CTIO_BVRIJHK(:,j), Trans_CTIO_BVRIJHK(:,j),&
                    Minmax_CTIO_BVRIJHK(:,j),size(Wavelength_CTIO_BVRIJHK(:,j)),ratio)
            flux = flux + ratio*simple_heating(i)*delta_mass(i)/(4*D**2)
        end do
        Magnitude_CTIO_BVRIJHK(j) = -2.5d0*log10(flux/ABzeropoint_CTIO_BVRIJHK(j))

    end do



    
    !CTIO_ugrizY
    do j = 1,6
        flux = 0.0d0
        if (photosphere_fell_on_the_center .eq. 0) then
            ! Blackbody assumption for the photosphere:
            call Blackbody_lambda_integral(T_eff,Wavelength_CTIO_ugrizY(:,j), Trans_CTIO_ugrizY(:,j),&
                    Minmax_CTIO_ugrizY(:,j),size(Wavelength_CTIO_ugrizY(:,j)),ratio)
            flux = flux + ratio*lum_photo/(4*D**2) ! not 4 pi D^2 because pi has been reduced
        end if
        ! Blackbody assumption for each radioactive layer above the photosphere:
        do i = index_photo,imax-1
            call Blackbody_lambda_integral(temp(i),Wavelength_CTIO_ugrizY(:,j), Trans_CTIO_ugrizY(:,j),&
                    Minmax_CTIO_ugrizY(:,j),size(Wavelength_CTIO_ugrizY(:,j)),ratio)
            flux = flux + ratio*simple_heating(i)*delta_mass(i)/(4*D**2)
        end do
        Magnitude_CTIO_ugrizY(j) = -2.5d0*log10(flux/ABzeropoint_CTIO_ugrizY(j))

    end do   

      open(666,file=trim(adjustl(outdir))//"/Magnitude_KNEC.dat",&
              status='unknown',position='append')
      write(666,"(13E18.9)") time, Magnitude_Gemini_ugriz,Magnitude_Gemini_JHKs
      write(666,"(13E18.9)") Magnitude_CTIO_BVRIJHK,Magnitude_CTIO_ugrizY
      close(666)

  endif


end subroutine analysis


! subroutine Filter(temp_x,Ratio_color,Ratio_color_mag,AB_zeropoint)
!     use blmod, only: Freq_UBVRI, Transmission_UBVRI
!     use parameters
!     use physical_constants
!     implicit none
!     real*8 :: temp_x
!     real*8 :: Ratio_color(5)
!     real*8 :: Ratio_color_mag(5)
!     real*8 :: AB_zeropoint(5)


!     real*8 :: Total_intensity
!     real*8 :: Color_intensity
!     real*8 :: Color_intensity_mag
!     integer :: i, j, n
!     real*8 :: b
!     real*8 :: df



!     Total_intensity = sigma_SB*temp_x**4/pi

!     do j = 1,5
!         Color_intensity = 0.0d0
!         n = size(Freq_UBVRI)
!         AB_zeropoint(j) = 0.0d0
!         ! integral intensity with filter
!         do i = 1,n
!             if (i .eq. 1) then
!                df = Freq_UBVRI(2)-Freq_UBVRI(1)
!                call Blackbody(temp_x,Freq_UBVRI(1),b)
!                Color_intensity = Color_intensity + 0.5d0*df*b*Transmission_UBVRI(1,j)
!                Color_intensity_mag = Color_intensity_mag + 0.5d0*df*b*Transmission_UBVRI(1,j)/(h_cgs*Freq_UBVRI(1))
!                AB_zeropoint(j) = AB_zeropoint(j) + 0.5d0*df*3631*Jy*Transmission_UBVRI(1,j)/(h_cgs*Freq_UBVRI(1))
!             elseif (i .eq. n) then
!                 df = Freq_UBVRI(n)-Freq_UBVRI(n-1)
!                 call Blackbody(temp_x,Freq_UBVRI(n),b)
!                 Color_intensity = Color_intensity + 0.5d0*df*b*Transmission_UBVRI(n,j)
!                 Color_intensity_mag = Color_intensity_mag + 0.5d0*df*b*Transmission_UBVRI(n,j)/(h_cgs*Freq_UBVRI(n))
!                 AB_zeropoint(j) = AB_zeropoint(j) + 0.5d0*df*3631*Jy*Transmission_UBVRI(n,j)/(h_cgs*Freq_UBVRI(n))

!             else
!                 df = Freq_UBVRI(i+1) - Freq_UBVRI(i-1)
!                 call Blackbody(temp_x,Freq_UBVRI(i),b)
!                 Color_intensity = Color_intensity + 0.5d0*df*b*Transmission_UBVRI(i,j)
!                 Color_intensity_mag = Color_intensity_mag + 0.5d0*df*b*Transmission_UBVRI(i,j)/(h_cgs*Freq_UBVRI(i))
!                 AB_zeropoint(j) = AB_zeropoint(j) + 0.5d0*df*3631*Jy*Transmission_UBVRI(i,j)/(h_cgs*Freq_UBVRI(i))

!             end if
!         end do
!         Ratio_color(j) = Color_intensity/Total_intensity
!         Ratio_color_mag(j) = Color_intensity_mag/Total_intensity
!     end do

! end subroutine Filter


! subroutine Blackbody(temp_x,nu,B)
!     use physical_constants
!     implicit none
!     real *8 :: temp_x
!     real *8 :: nu
!     real *8 :: B
!     B = 2*h_cgs*nu**3/clite**2*(exp(h_cgs*nu/kboltz/temp_x)-1)**(-1)

! end subroutine Blackbody

subroutine Blackbody_lambda(temp_x,lambda,B)
!   2c/(lambda^4 sigma T^4)/(exp(hc/(lambda k T))-1), for magnitude calculation
    use physical_constants
    implicit none
    real *8 :: temp_x
    real *8 :: lambda
    real *8 :: B
    B = 2*clite/lambda**4/sigma_SB/temp_x**4*(exp(h_cgs*clite/lambda/kboltz/temp_x)-1)**(-1)

end subroutine Blackbody_lambda

subroutine Blackbody_lambda_integral(temp_x,filter_lambda,filter_transmission,filter_minmax,filter_size,I)
    !integrate Blackbody_lambda with filter
    use physical_constants
    implicit none
    real*8 :: temp_x
    integer :: filter_size
    real*8 :: filter_lambda(filter_size),filter_transmission(filter_size)
    real*8 :: filter_minmax(2) !min(wavelength)=filter_minmax(1),max(wavelength)=filter_minmax(2)

    !output:
    real*8 :: I

    !local:
    integer :: n_interval = 100 !n>2 and n should be even number
    real*8 :: h, B, t
    integer :: j
    ! Use Composite Simpson's rule

    I = 0.0d0
    h = (filter_minmax(2)-filter_minmax(1))/n_interval
    call Blackbody_lambda(temp_x,filter_minmax(1),B)
    call map_map(t,filter_minmax(1),filter_transmission,filter_lambda,filter_size)
    I = I + B*t

    call Blackbody_lambda(temp_x,filter_minmax(2),B)
    call map_map(t,filter_minmax(2),filter_transmission,filter_lambda,filter_size)
    I = I + B*t


    do j = 1,n_interval/2
        call Blackbody_lambda(temp_x,filter_minmax(1)+(2*j-1)*h,B)
        call map_map(t,filter_minmax(1)+(2*j-1)*h,filter_transmission,filter_lambda,filter_size)
        I = I + 4*B*t
    end do
    do j = 1,n_interval/2-1
        call Blackbody_lambda(temp_x,filter_minmax(1)+(2*j)*h,B)
        call map_map(t,filter_minmax(1)+(2*j)*h,filter_transmission,filter_lambda,filter_size)
        I = I + 2*B*t
    end do
    I = I*h/3

end subroutine Blackbody_lambda_integral


