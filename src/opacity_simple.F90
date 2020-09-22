subroutine opacity_simple(temp_x,kappa_x,kappa_table_x,dkappadt_x)

    use blmod, only: ye,logT
    use parameters
    use physical_constants
    implicit none

    !input:
    real*8 temp_x(imax)

    !output:
    real*8 kappa_x(imax)
    real*8 kappa_table_x(imax)
    real*8 dkappadt_x(imax)

    !local:
    integer i
    real*8 :: kappa_min = 1.0d0
    real*8 :: kappa_max = 10.0d0


    do i=1, imax - 1
        logT(i)=log10(temp_x(i))

        if (ye(i) .gt. 0.25) then
            kappa_x(i)= kappa_min
        else
            kappa_x(i) = kappa_max
        end if


        !derivative is not used in the current version of the code
        dkappadt_x(i) = 0.0d0

    end do

    kappa_table_x(1:imax-1) = kappa_x(1:imax-1)


    ! Since the temperature temp(imax) is evaluated in the evolution, but
    ! set by the boundary condition, it doesn't make sense to find
    ! kappa(imax) from the table, so we assume:
    kappa_x(imax) = kappa_x(imax-1)
    kappa_table_x(imax) = kappa_table_x(imax-1)
    dkappadt_x(imax) = dkappadt_x(imax-1)

end subroutine opacity_simple
