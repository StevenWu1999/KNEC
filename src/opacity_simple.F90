subroutine opacity_simple(temp_x,kappa_x,kappa_table_x,dkappadt_x)

    use blmod, only: ye_initial,logT
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
        kappa_x(i) = 1.0d0 + 9.0d0/(1.0d0+(4.0d0*ye_initial(i))**5)

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

