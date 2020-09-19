subroutine nuclear_heating_rate
    use blmod, only : ye, heating_rate_fitting, heating_deposit_function, simple_heating,&
            r, rho, time, delta_mass
    use parameters
    use physical_constants
    implicit none

    !local:
    !Ye correction
    integer :: i
    real*8 :: ye_correction(imax)
    real*8, parameter :: epsilon_min = 0.5
    real*8, parameter :: epsilon_max = 2.5
    real*8, parameter :: t_eps = 1 * 24.0 * 3600.0 !1day
    !radiation transfer
    real*8 :: r_heating
    real*8 :: delta_r
    real*8 :: r_x, r_max
    real*8 :: r_j
    real*8 :: th, delta_th
    real*8 :: th_min(imax)
    real*8, parameter :: th_max = pi
    real*8 :: delta_tau_j
    real*8 :: I_prime, I_prime_av
    integer, parameter :: npoints_radial_integration = 150
    integer, parameter :: npoints_angular_integration = 150
    integer :: index, ibuffer

    ! correction for Ye>0.25:

    do i = 1, imax
        if (ye(i) .gt. 0.25) then
            ye_correction(i) = epsilon_min + epsilon_max * (1 + exp(4 * (time / t_eps - 1)))**(-1)
        else
            ye_correction(i) = 1.0
        end if
    end do

    heating_rate_fitting = heating_epsilon_0 * ye_correction * (heating_epsilon_th / 0.5) &
            * (0.5 - 1.0 / pi * atan((time - heating_t_0) / heating_sigma))

    ! radiation transfer
    r_heating = r(imax)

    th_min(1:imax - 1) = 0.0d0
    th_min(imax) = pi * 0.5d0

    do i = 1, imax
        th = th_max
        I_prime_av = 0
        delta_th = (th_max - th_min(i)) / npoints_angular_integration

        do while(th.gt.(th_min(i) + 1.0d-14))
            r_max = -r(i) * cos(th) + sqrt((r(i) * cos(th))**2 - (r(i)**2 - r_heating**2))
            delta_r = r_max / npoints_radial_integration
            I_prime = 0
            r_x = r_max

            do while(r_x.gt.0)
                r_j = sqrt(r(i) * r(i) + r_x * r_x + 2.0d0 * r(i) * r_x * cos(th))

                if(r_j.le.r(1)) then !inside the excised region
                    delta_tau_j = 0.0d0
                else if(r_j.ge.r(imax - 1)) then
                    delta_tau_j = delta_r * ye(imax - 1) * 0.06d0 * rho(imax - 1)
                else
                    call map_find_index(imax, r, r_j, ibuffer, index)
                    delta_tau_j = delta_r * ye(index) * 0.06d0 * rho(index)
                end if

                I_prime = (I_prime - 1) * exp(-delta_tau_j) + 1
                r_x = r_x - delta_r
            end do

            I_prime_av = I_prime_av + I_prime * sin(th) * delta_th * 0.5d0
            th = th - delta_th
        end do

        heating_deposit_function(i) = I_prime_av
    end do

    simple_heating(1:imax)=heating_deposit_function(1:imax)*heating_rate_fitting(1:imax)

end subroutine nuclear_heating_rate

