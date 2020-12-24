! Created by Wu Zhenyu on 11/22/20.

subroutine read_filter
    use blmod, only: Freq_UBVRI, Transmission_UBVRI,&
            Wavelength_Gemini_ugriz,Trans_Gemini_ugriz,Minmax_Gemini_ugriz,ABzeropoint_Gemini_ugriz,&
            Wavelength_Gemini_JHKs,Trans_Gemini_JHKs,Minmax_Gemini_JHKs,ABzeropoint_Gemini_JHKs,&
            Wavelength_CTIO_BVRIJHK,Trans_CTIO_BVRIJHK,Minmax_CTIO_BVRIJHK,ABzeropoint_CTIO_BVRIJHK,&
            Wavelength_CTIO_ugrizY,Trans_CTIO_ugrizY,Minmax_CTIO_ugrizY,ABzeropoint_CTIO_ugrizY
    use parameters
    implicit none

    integer :: nlines, i , j
    real*8 :: test,zeropoint
    character*200 :: linestring

    open(666,file=trim("filter/Johnson-Cousins_UBVRI/UBVRI.dat"),status='unknown',&
            form='formatted',action='read')

    read(666,*)
    read(666,*) nlines

    allocate(Freq_UBVRI(nlines))
    allocate(Transmission_UBVRI(nlines,5))
    do i = 1,nlines
        read(666,*) Freq_UBVRI(i), Transmission_UBVRI(i,1:5)
    end do
    close(666)


    ! Gemini GMOS-N ugriz
    open(666,file=trim("filter/Gemini/Gemini_GMOS-N_ugriz.dat"),status='unknown',&
            form='formatted',action='read')
    do i = 1,3
        read(666,*)
    end do
    read(666,*) Minmax_Gemini_ugriz(1,:)
    read(666,*) Minmax_Gemini_ugriz(2,:)
    read(666,*) 
    read(666,*) nlines
    allocate(Wavelength_Gemini_ugriz(nlines,5))
    allocate(Trans_Gemini_ugriz(nlines,5))
    do i = 1,nlines
        read(666,*) Wavelength_Gemini_ugriz(i,1), Trans_Gemini_ugriz(i,1),&
                Wavelength_Gemini_ugriz(i,2), Trans_Gemini_ugriz(i,2),&
                Wavelength_Gemini_ugriz(i,3), Trans_Gemini_ugriz(i,3),&
                Wavelength_Gemini_ugriz(i,4), Trans_Gemini_ugriz(i,4),&
                Wavelength_Gemini_ugriz(i,5), Trans_Gemini_ugriz(i,5)
    end do
    close(666)

    do i=1,5
        call Band_ABzeropoint(Trans_Gemini_ugriz(:,i),Wavelength_Gemini_ugriz(:,i),&
                Minmax_Gemini_ugriz(:,i),size(Wavelength_Gemini_ugriz(:,i)),ABzeropoint_Gemini_ugriz(i))
    end do



   ! Gemini Flamingos2 JHKs
    open(666,file=trim("filter/Gemini/Gemini_Flamingos2_JHKs.dat"),status='unknown',&
            form='formatted',action='read')
    do i = 1,3
        read(666,*)
    end do
    read(666,*) Minmax_Gemini_JHKs(1,:)
    read(666,*) Minmax_Gemini_JHKs(2,:)
    read(666,*)
    read(666,*) nlines
    allocate(Wavelength_Gemini_JHKs(nlines,3))
    allocate(Trans_Gemini_JHKs(nlines,3))
    do i = 1,nlines
        read(666,*) Wavelength_Gemini_JHKs(i,1), Trans_Gemini_JHKs(i,1),&
                Wavelength_Gemini_JHKs(i,2), Trans_Gemini_JHKs(i,2),&
                Wavelength_Gemini_JHKs(i,3), Trans_Gemini_JHKs(i,3)
    end do
    close(666)

    do i=1,3
        call Band_ABzeropoint(Trans_Gemini_JHKs(:,i),Wavelength_Gemini_JHKs(:,i),&
                Minmax_Gemini_JHKs(:,i),size(Wavelength_Gemini_JHKs(:,i)),ABzeropoint_Gemini_JHKs(i))
    end do


    ! CTIO Andicam BVRIJHK
    open(666,file=trim("filter/CTIO/CTIO_Andicam_BVRIJHK.dat"),status='unknown',&
            form='formatted',action='read')
    do i = 1,3
        read(666,*)
    end do
    read(666,*) Minmax_CTIO_BVRIJHK(1,:)
    read(666,*) Minmax_CTIO_BVRIJHK(2,:)
    read(666,*)
    read(666,*) nlines
    allocate(Wavelength_CTIO_BVRIJHK(nlines,7))
    allocate(Trans_CTIO_BVRIJHK(nlines,7))
    do i = 1,nlines
        read(666,*) Wavelength_CTIO_BVRIJHK(i,1), Trans_CTIO_BVRIJHK(i,1),&
                Wavelength_CTIO_BVRIJHK(i,2), Trans_CTIO_BVRIJHK(i,2),&
                Wavelength_CTIO_BVRIJHK(i,3), Trans_CTIO_BVRIJHK(i,3),&
                Wavelength_CTIO_BVRIJHK(i,4), Trans_CTIO_BVRIJHK(i,4),&
                Wavelength_CTIO_BVRIJHK(i,5), Trans_CTIO_BVRIJHK(i,5),&
                Wavelength_CTIO_BVRIJHK(i,6), Trans_CTIO_BVRIJHK(i,6),&
                Wavelength_CTIO_BVRIJHK(i,7), Trans_CTIO_BVRIJHK(i,7)
    end do
    close(666)

    do i=1,7
        call Band_ABzeropoint(Trans_CTIO_BVRIJHK(:,i),Wavelength_CTIO_BVRIJHK(:,i),&
                Minmax_CTIO_BVRIJHK(:,i),size(Wavelength_CTIO_BVRIJHK(:,i)),ABzeropoint_CTIO_BVRIJHK(i))
    end do


    ! CTIO DECam ugrizY
    open(666,file=trim("filter/CTIO/CTIO_DECam_ugrizY.dat"),status='unknown',&
            form='formatted',action='read')
    do i = 1,3
        read(666,*)
    end do
    read(666,*) Minmax_CTIO_ugrizY(1,:)
    read(666,*) Minmax_CTIO_ugrizY(2,:)
    read(666,*)
    read(666,*) nlines
    allocate(Wavelength_CTIO_ugrizY(nlines,6))
    allocate(Trans_CTIO_ugrizY(nlines,6))
    do i = 1,nlines
        read(666,*) Wavelength_CTIO_ugrizY(i,1), Trans_CTIO_ugrizY(i,1),&
                Wavelength_CTIO_ugrizY(i,2), Trans_CTIO_ugrizY(i,2),&
                Wavelength_CTIO_ugrizY(i,3), Trans_CTIO_ugrizY(i,3),&
                Wavelength_CTIO_ugrizY(i,4), Trans_CTIO_ugrizY(i,4),&
                Wavelength_CTIO_ugrizY(i,5), Trans_CTIO_ugrizY(i,5),&
                Wavelength_CTIO_ugrizY(i,6), Trans_CTIO_ugrizY(i,6)
    end do
    close(666)

    do i=1,6
        call Band_ABzeropoint(Trans_CTIO_ugrizY(:,i),Wavelength_CTIO_ugrizY(:,i),&
                Minmax_CTIO_ugrizY(:,i),size(Wavelength_CTIO_ugrizY(:,i)),ABzeropoint_CTIO_ugrizY(i))
    end do



    open(666,file=trim(adjustl(outdir))//"/Magnitude_KNEC.dat",&
            status='unknown',position='append')
    linestring="time[s]  Gemini_u_g_r_i_z  Gemini_J_H_Ks"
    write(666,*) linestring
    linestring="CTIO_B_V_R_I_J_H_K  CTIO_u_g_r_i_z_Y"
    write(666,*) linestring
    close(666)




end subroutine read_filter




subroutine Band_ABzeropoint(filter_transmission,filter_lambda,filter_minmax,filter_size,zeropoint)
    use physical_constants
    implicit none
    integer :: filter_size
    real*8 :: filter_lambda(filter_size),filter_transmission(filter_size)
    real*8 :: filter_minmax(2) !min(wavelength)=filter_minmax(1),max(wavelength)=filter_minmax(2)

    integer :: n_interval = 10000 !n>2 and n should be even number
    real*8 :: h, t, lambda
    integer :: j

    !output:
    real*8 :: zeropoint
    zeropoint = 0.0d0
    ! Composite Simpson's rule
    h = (filter_minmax(2)-filter_minmax(1))/n_interval
    call map_map(t,filter_minmax(1),filter_transmission,filter_lambda,filter_size)
    zeropoint = zeropoint + 3631*Jy/h_cgs/filter_minmax(1)*t

    call map_map(t,filter_minmax(2),filter_transmission,filter_lambda,filter_size)
    zeropoint = zeropoint + 3631*Jy/h_cgs/filter_minmax(2)*t

    do j = 1,n_interval/2
        lambda = filter_minmax(1)+(2*j-1)*h
        call map_map(t,lambda,filter_transmission,filter_lambda,filter_size)
        zeropoint = zeropoint + 4*3631*Jy/h_cgs/lambda*t
    end do

    do j = 1,n_interval/2-1
        lambda = filter_minmax(1)+(2*j)*h
        call map_map(t,lambda,filter_transmission,filter_lambda,filter_size)
        zeropoint = zeropoint + 2*3631*Jy/h_cgs/lambda*t
    end do
    zeropoint = zeropoint*h/3


end subroutine Band_ABzeropoint