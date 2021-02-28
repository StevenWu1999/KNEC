subroutine nuclear_heating_rate
    use blmod, only : ye_initial ,entropy_frominput, heating_rate_fitting, heating_deposit_function, simple_heating,&
            r, rho, time, delta_mass, A_int_ricigliano, expansion_timescale, A_int_arctan, B_int_arctan
    use parameters
    use physical_constants
    use heating_rate_LR15_module
    use heating_rate_arctan_module
    implicit none

    !local:
    !Ye correction
    integer :: i
    real*8 :: ye_correction(imax)
    real*8, parameter :: epsilon_min = 0.5
    real*8, parameter :: epsilon_max = 2.5
    real*8, parameter :: t_eps = 1 * 24.0 * 3600.0 !1day
    real*8 :: t_day
    real*8 :: hybrid_K, hybrid_R, K_fraction, R_fraction


    if (trim(adjustl(heating_formula)) .eq. "Korobkin") then

        ! correction for Ye>0.25:
        do i = 1, imax
            if (ye_initial(i) .gt. 0.25) then
                ye_correction(i) = epsilon_min + epsilon_max * (1 + exp(4 * (time / t_eps - 1)))**(-1)
            else
                ye_correction(i) = 1.0d0
            end if
        end do

        heating_rate_fitting = heating_epsilon_0 * ye_correction * (heating_epsilon_th / 0.5) &
                * (0.5 - 1.0 / pi * atan((time - heating_t_0) / heating_sigma))**heating_alpha
        simple_heating(:) = heating_rate_fitting(:)

    else if (trim(adjustl(heating_formula)) .eq. "LR15") then

        do i = 1,imax
            call calc_heating_rate_LR15(expansion_timescale(i),entropy_frominput(i),ye_initial(i),time,simple_heating(i))
        end do
        simple_heating(:) = heating_epsilon_th*simple_heating(:)

    else if (trim(adjustl(heating_formula)) .eq. "Ricigliano") then
        t_day = time/(24.0*3600.0)
        do i = 1,imax
            simple_heating(i) = heating_epsilon_th*A_int_ricigliano(i,1)*t_day**(-A_int_ricigliano(i,2))
        end do

    else if (trim(adjustl(heating_formula)) .eq. "Hybrid") then

        if (time .le. 0.9d0*heating_t_cutoff) then
            !Korobkin
            ! correction for Ye>0.25:
            do i = 1, imax
                if (ye_initial(i) .gt. 0.25) then
                    ye_correction(i) = epsilon_min + epsilon_max * (1 + exp(4 * (time / t_eps - 1)))**(-1)
                else
                    ye_correction(i) = 1.0d0
                end if
            end do

            simple_heating(:) = heating_epsilon_0 * ye_correction * (heating_epsilon_th / 0.5) &
                    * (0.5 - 1.0 / pi * atan((time - heating_t_0) / heating_sigma))**heating_alpha

        else if (time .ge. 1.1d0*heating_t_cutoff) then
            ! Ricigliano
            t_day = time/(24.0*3600.0)
            do i = 1,imax
                simple_heating(i) = heating_epsilon_th*A_int_ricigliano(i,1)*t_day**(-A_int_ricigliano(i,2))
            end do
        else
            do i = 1, imax
                if (ye_initial(i) .gt. 0.25) then
                    ye_correction(i) = epsilon_min + epsilon_max * (1 + exp(4 * (time / t_eps - 1)))**(-1)
                else
                    ye_correction(i) = 1.0d0
                end if
                hybrid_K = heating_epsilon_0 * ye_correction(i) * (heating_epsilon_th / 0.5) &
                        * (0.5 - 1.0 / pi * atan((time - heating_t_0) / heating_sigma))**heating_alpha

                t_day = time/(24.0*3600.0)
                hybrid_R = heating_epsilon_th*A_int_ricigliano(i,1)*t_day**(-A_int_ricigliano(i,2))


                R_fraction = (time-0.9d0*heating_t_cutoff)/(0.2d0*heating_t_cutoff)

                K_fraction = 1.0d0 - R_fraction
                simple_heating(i) = K_fraction*hybrid_K + R_fraction * hybrid_R

            end do



        end if

    elseif (trim(adjustl(heating_formula)) .eq. "arctan") then
        t_day = time/(24.0*3600.0)
        do i = 1,imax

            if (t_day.ge.arctan_t2) then
                simple_heating(i) = A_int_arctan(i,1)*t_day**(-A_int_arctan(i,2))
            else
                simple_heating(i) = B_int_arctan(i,1)*(0.5-oneoverpi*datan((time-arctan_t0)/arctan_sigma0))**B_int_arctan(i,2)
            end if
        end do

    end if







end subroutine nuclear_heating_rate

