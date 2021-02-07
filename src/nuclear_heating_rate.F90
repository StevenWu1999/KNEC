subroutine nuclear_heating_rate
    use blmod, only : ye_initial ,entropy_frominput, heating_rate_fitting, heating_deposit_function, simple_heating,&
            r, rho, time, delta_mass, A_int_ricigliano, expansion_timescale
    use parameters
    use physical_constants
    use heating_rate_LR15_module
    implicit none

    !local:
    !Ye correction
    integer :: i
    real*8 :: ye_correction(imax)
    real*8, parameter :: epsilon_min = 0.5
    real*8, parameter :: epsilon_max = 2.5
    real*8, parameter :: t_eps = 1 * 24.0 * 3600.0 !1day
    real*8 :: t_day


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

    end if






end subroutine nuclear_heating_rate

