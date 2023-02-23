module helmholtz

  ! here is the tabular helmholtz free energy eos:
  !
  ! routine read_helm_table reads an electron helm free energy table
  ! routine helmholtz_eos computes the pressure, energy and entropy via tables

  implicit none

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%% Declare variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! sizes of the tables
  ! normal table, big table, bigger table, denser bigger table
  integer ::       imax_helm,jmax_helm
  !parameter        (imax_helm = 211, jmax_helm = 71)    ! original
  parameter        (imax_helm = 271, jmax_helm = 101)    ! standard
  ! parameter        (imax_helm = 541, jmax_helm = 201)  ! twice as dense
  ! parameter        (imax_helm = 136, jmax_helm = 51)   ! half as dense

  ! for the electrons table
  ! density and temperature
  double precision :: d_table(imax_helm),t_table(jmax_helm)
  double precision :: tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi
  ! for standard table limits
  parameter           (tlo   = 3.0d0, &
                       thi   = 13.0d0, &
                       tstp  = (thi - tlo)/float(jmax_helm-1), &
                       tstpi = 1.0d0/tstp, &
                       dlo   = -12.0d0, &
                       dhi   = 15.0d0, &
                       dstp  = (dhi - dlo)/float(imax_helm-1), &
                       dstpi = 1.0d0/dstp)
  ! for the helmholtz free energy tables
  double precision :: f(imax_helm,jmax_helm),fd(imax_helm,jmax_helm), &
                      ft(imax_helm,jmax_helm),fdd(imax_helm,jmax_helm), &
                      ftt(imax_helm,jmax_helm), fdt(imax_helm,jmax_helm), &
                      fddt(imax_helm,jmax_helm),fdtt(imax_helm,jmax_helm), &
                      fddtt(imax_helm,jmax_helm)
  ! for the pressure derivative with density tables
  double precision :: dpdf(imax_helm,jmax_helm),dpdfd(imax_helm,jmax_helm), &
                      dpdft(imax_helm,jmax_helm),dpdfdt(imax_helm,jmax_helm)
  ! for chemical potential tables
  double precision :: ef(imax_helm,jmax_helm),efd(imax_helm,jmax_helm), &
                      eft(imax_helm,jmax_helm),efdt(imax_helm,jmax_helm)
  ! for the number density tables
  double precision :: xf(imax_helm,jmax_helm),xfd(imax_helm,jmax_helm), &
                      xft(imax_helm,jmax_helm),xfdt(imax_helm,jmax_helm)
  ! for storing the differences
  double precision :: dt_sav(jmax_helm),dt2_sav(jmax_helm), &
                      dti_sav(jmax_helm),dt2i_sav(jmax_helm), &
                      dt3i_sav(jmax_helm), dd_sav(imax_helm), &
                      dd2_sav(imax_helm), ddi_sav(imax_helm), &
                      dd2i_sav(imax_helm),dd3i_sav(imax_helm)

  ! helhmoltz table path
  character(80) :: helm_table_path
  parameter        (helm_table_path = "./tables/helm_table.dat")
  ! check if the helhmoltz table has been read
  logical :: if_table_read = .false.
  ! failure of the eos
  logical :: eosfail = .false.

  contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%% Routine to read Helmholtz table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine read_helm_table()

      ! this routine reads the helmholtz eos file, and
      ! MUST be called once before the helmholtz_eos routine is invoked

      implicit none
      ! save

      ! declare local variables
      integer :: i,j
      double precision :: tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                          dd,dd2,ddi,dd2i,dd3i

      ! open the file (use softlinks to input the desired table)
      open(unit=19, file=trim(helm_table_path), status='old')

      !.....Read Helmholtz EoS tabulated data....................................
      ! read the helmholtz free energy and its derivatives (this information
      ! is stored in the first block of rows from the helm_table file)
      do j=1,jmax_helm
        tsav = tlo + (j-1)*tstp
        t_table(j) = 10.0d0**(tsav)
        do i=1,imax_helm
            dsav = dlo + (i-1)*dstp
            d_table(i) = 10.0d0**(dsav)
            read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                       fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
      enddo
      ! read the pressure derivative with density table
      do j=1,jmax_helm
         do i=1,imax_helm
            read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
         enddo
      enddo
      ! read the electron chemical potential table
      do j=1,jmax_helm
        do i=1,imax_helm
            read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
      enddo
      ! read the number density table
      do j=1,jmax_helm
        do i=1,imax_helm
            read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
      enddo
      ! close the file and write a summary message
      close(unit=19)

      ! Construct the temperature and density deltas and their inverses
      do j=1,jmax_helm-1
        dth         = t_table(j+1) - t_table(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt3i        = dt2i*dti
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
        dt3i_sav(j) = dt3i
      end do
      do i=1,imax_helm-1
        dd          = d_table(i+1) - d_table(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd3i        = dd2i*ddi
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
        dd3i_sav(i) = dd3i
      enddo

      if_table_read = .true.

    end subroutine read_helm_table


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%% Routine to compute thermodynamic quantities %%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine helmholtz_eos(rhox,tempx,yex,abarx,px,ex,sx,cs2x,dpdtx,dedtx,pradx,np)

      ! given a density den [g/cm**3], temperature temp [K], and a composition
      ! characterized by abar and zbar, this routine returns most of the other
      ! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
      ! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
      ! their derivatives with respect to temperature, density, abar, and zbar.
      ! other quantites such the normalized chemical potential eta (plus its
      ! derivatives), number density of electrons and positron pair (along
      ! with their derivatives), adiabatic indices, specific heats, and
      ! relativistically correct sound speed are also returned.
      !
      ! this routine assumes planckian photons, an ideal gas of ions,
      ! and an electron-positron gas with an arbitrary degree of relativity
      ! and degeneracy. interpolation in a table of the helmholtz free energy
      ! is used to return the electron-positron thermodynamic quantities.
      ! all other derivatives are analytic.
      !
      ! references: cox & giuli chapter 24 ; timmes & swesty apj 1999
      !
      ! subroutine read_helm_table MUST has been called once before calling this
      ! eos subroutine

      use physical_constants
      implicit none
      ! save

      !.....Declare..............................................................

      ! input:
      integer :: np
      real*8 :: rhox(np),tempx(np),yex(np),abarx(np)
      ! output:
      real*8 :: px(np),ex(np),sx(np),cs2x(np),dpdtx(np),dedtx(np),pradx(np)

      integer :: i,j
      double precision :: temp,den,abar,zbar,ytot1,ye, &
                          x,y,zz,zzi,deni,tempi,xni,dxnidd, &
                          dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                          dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                          deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                          dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                          sion,pele,eele,sele,pres,ener,entr,dpresdd, & !,xnem
                          dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv, &
                          gam1,chit,chid,sound,s
      double precision :: pgas,dpgasdd,dpgasdt, &
                          egas,degasdd,degasdt, &
                          sgas,dsgasdd,dsgasdt

      double precision :: sioncon,forth,kergavo,asoli3,light2
      parameter           (sioncon = (2.0d0 * pi * amu * kboltz)/(h_cgs*h_cgs), &
                           forth   = 4.0d0/3.0d0, &
                           kergavo = kboltz * avo_real, &
                           asoli3  = a_rad/3.0d0, &
                           light2  = clite * clite)

      ! for the interpolations
      integer          :: iat,jat
      double precision :: free,df_d,df_t,df_tt,df_dt!,df_dd
      double precision :: xt,xd,mxt,mxd, &
                          si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                          si0d,si1d,si2d,si0md,si1md,si2md, &
                          dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                          dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                          ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                          z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
                          dpsi2,ddpsi2,din,h5,fi(36), &
                          xpsi0,xpsi1,h3, &
                          w0t,w1t,w2t,w0mt,w1mt,w2mt, &
                          w0d,w1d,w2d,w0md,w1md,w2md

      ! for the uniform background coulomb correction
      double precision :: dsdd,lami,inv_lami,lamidd, &
                          plasg,plasgdd,plasgdt, &
                          ecoul,decouldd,decouldt, &
                          pcoul,dpcouldd,dpcouldt, &
                          scoul,dscouldd,dscouldt, &
                          a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
      parameter           (a1    = -0.898004d0, &
                           b1    =  0.96786d0, &
                           c1    =  0.220703d0, &
                           d1    = -0.86097d0, &
                           e1    =  2.5269d0, &
                           a2    =  0.29561d0, &
                           b2    =  1.9885d0, &
                           c2    =  0.288675d0, &
                           third =  1.0d0/3.0d0, &
                           esqu  =  qe * qe)

      ! quintic hermite polynomial statement functions
      ! psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)

      ! psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)

      ! psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)

      ! biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
             fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt
      ! cubic hermite polynomial statement functions
      ! psi0
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      ! psi1
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      ! bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
             fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
           + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
           + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
           + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
           + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt

      ! popular format statement
      01    format(1x,5(a,1pe11.3))

      !..........................................................................
      !..........................................................................
      !.....Start of pipeline loop, normal execution starts here.................
      !..........................................................................
      !..........................................................................

      if (if_table_read .eqv. .false.) then
        write(*,*) "Helmholtz eos called before reading Helhmoltz table! exit!"
        call exit(-1)
      end if

      eosfail = .false.
      do j=1,np
        ! if (tempx(j) .le. 0.0) stop 'temp less than 0 in helmeos'
        ! if (rhox(j)  .le. 0.0) stop 'den less than 0 in helmeos'
        temp  = tempx(j)
        den   = rhox(j)
        ye = yex(j)
        abar  = abarx(j)
        zbar = abar*ye
        ytot1 = 1.0d0/abar
        ! Initialize
        deni    = 1.0d0/den
        tempi   = 1.0d0/temp
        kt      = kboltz * temp
        ktinv   = 1.0d0/kt

        !......................................................................
        !.....Radiation section: planckian photons.............................
        !......................................................................

        ! pressure & derivatives
        prad    = asoli3 * temp * temp * temp * temp
        dpraddd = 0.0d0
        dpraddt = 4.0d0 * prad*tempi

        ! specific thermal energy & derivatives
        erad    = 3.0d0 * prad*deni
        deraddd = -erad*deni
        deraddt = 3.0d0 * dpraddt*deni

        ! entropy & derivatives
        srad    = (prad*deni + erad)*tempi
        dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
        dsraddt = (dpraddt*deni + deraddt - srad)*tempi

        !......................................................................
        !.....Ion section: ideal gas...........................................
        !......................................................................

        ! ion number density & derivatives
        xni     = avo_real * ytot1 * den
        dxnidd  = avo_real * ytot1

        ! pressure & derivatives
        pion    = xni * kt
        dpiondd = dxnidd * kt
        dpiondt = xni * kboltz

        ! specific thermal energy & derivatives
        eion    = 1.5d0 * pion*deni
        deiondd = (1.5d0 * dpiondd - eion)*deni
        deiondt = 1.5d0 * dpiondt*deni

        ! sackur-tetrode equation for the ion entropy of
        ! a single ideal gas characterized by abar
        ! the implemented expression generalizes Lippuner & Roberts 2018, eq (A15)
        ! for a monoatomic ideal gas, G(T)=1 and (pion*deni+eion)*tempi*m =5/2*kboltz
        x = abar*abar*sqrt(abar) * deni/avo_real
        s = sioncon * temp
        z = x * s * sqrt(s)
        y = log(z)
        ! entropy & derivatives
        sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
        dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                   - kergavo * deni * ytot1
        dsiondt = (dpiondt*deni + deiondt)*tempi - &
                  (pion*deni + eion) * tempi*tempi &
                  + 1.5d0 * kergavo * tempi*ytot1

        !......................................................................
        !.....Electron-positron section: helmholtz free energy table...........
        !......................................................................

        ! electron number density: assume complete ionization: not used in kNEC
        ! xnem = xni * zbar

        !..................................................................
        !.....Read Helmholtz table.........................................

        ! enter the table with the folowing value for the proton mass density din
        ! (remember that ye = zbar/abar -> ye*den = n*zbar*amu):
        din = ye*den

        ! check input boundaries
        if (temp .gt. t_table(jmax_helm)) then
            write(6,01) '   temp=',temp,' t_table(jmax_helm)=',t_table(jmax_helm)
!            write(6,*) 'temp too hot, off grid'
!            write(6,*) 'setting eosfail to true and returning'
            eosfail = .true.
            return
        end if
        if (temp .lt. t_table(1)) then
            write(6,01) '   temp=',temp,' t_table(1)=',t_table(1)
!            write(6,*) 'temp too cold, off grid'
!            write(6,*) 'setting eosfail to true and returning'
            eosfail = .true.
            return
        end if
        if (din  .gt. d_table(imax_helm)) then
            write(6,01) '   den*ye=',din,' d_table(imax_helm)=',d_table(imax_helm)
!            write(6,*) 'ye*den too big, off grid'
!            write(6,*) 'setting eosfail to true and returning'
            eosfail = .true.
            return
        end if
        if (din  .lt. d_table(1)) then
            write(6,01) '   ye*den=',din,' d_table(1)=',d_table(1)
!            write(6,*) 'ye*den too small, off grid'
!            write(6,*) 'setting eosfail to true and returning'
            eosfail = .true.
            return
        end if

        !.....Compute free energy and derivatives................

        ! hash locate this temperature and density
        jat = int((log10(temp) - tlo)*tstpi) + 1
        jat = max(1,min(jat,jmax_helm-1))
        iat = int((log10(din) - dlo)*dstpi) + 1
        iat = max(1,min(iat,imax_helm-1))

        ! access the table locations only once
        fi(1)  = f(iat,jat)
        fi(2)  = f(iat+1,jat)
        fi(3)  = f(iat,jat+1)
        fi(4)  = f(iat+1,jat+1)
        fi(5)  = ft(iat,jat)
        fi(6)  = ft(iat+1,jat)
        fi(7)  = ft(iat,jat+1)
        fi(8)  = ft(iat+1,jat+1)
        fi(9)  = ftt(iat,jat)
        fi(10) = ftt(iat+1,jat)
        fi(11) = ftt(iat,jat+1)
        fi(12) = ftt(iat+1,jat+1)
        fi(13) = fd(iat,jat)
        fi(14) = fd(iat+1,jat)
        fi(15) = fd(iat,jat+1)
        fi(16) = fd(iat+1,jat+1)
        fi(17) = fdd(iat,jat)
        fi(18) = fdd(iat+1,jat)
        fi(19) = fdd(iat,jat+1)
        fi(20) = fdd(iat+1,jat+1)
        fi(21) = fdt(iat,jat)
        fi(22) = fdt(iat+1,jat)
        fi(23) = fdt(iat,jat+1)
        fi(24) = fdt(iat+1,jat+1)
        fi(25) = fddt(iat,jat)
        fi(26) = fddt(iat+1,jat)
        fi(27) = fddt(iat,jat+1)
        fi(28) = fddt(iat+1,jat+1)
        fi(29) = fdtt(iat,jat)
        fi(30) = fdtt(iat+1,jat)
        fi(31) = fdtt(iat,jat+1)
        fi(32) = fdtt(iat+1,jat+1)
        fi(33) = fddtt(iat,jat)
        fi(34) = fddtt(iat+1,jat)
        fi(35) = fddtt(iat,jat+1)
        fi(36) = fddtt(iat+1,jat+1)

        ! various differences
        xt  = max( (temp - t_table(jat))*dti_sav(jat), 0.0d0)
        xd  = max( (din - d_table(iat))*ddi_sav(iat), 0.0d0)
        mxt = 1.0d0 - xt
        mxd = 1.0d0 - xd

        ! the six density and six temperature basis functions
        si0t =   psi0(xt)
        si1t =   psi1(xt)*dt_sav(jat)
        si2t =   psi2(xt)*dt2_sav(jat)

        si0mt =  psi0(mxt)
        si1mt = -psi1(mxt)*dt_sav(jat)
        si2mt =  psi2(mxt)*dt2_sav(jat)

        si0d =   psi0(xd)
        si1d =   psi1(xd)*dd_sav(iat)
        si2d =   psi2(xd)*dd2_sav(iat)

        si0md =  psi0(mxd)
        si1md = -psi1(mxd)*dd_sav(iat)
        si2md =  psi2(mxd)*dd2_sav(iat)

        ! derivatives of the weight functions
        dsi0t =   dpsi0(xt)*dti_sav(jat)
        dsi1t =   dpsi1(xt)
        dsi2t =   dpsi2(xt)*dt_sav(jat)

        dsi0mt = -dpsi0(mxt)*dti_sav(jat)
        dsi1mt =  dpsi1(mxt)
        dsi2mt = -dpsi2(mxt)*dt_sav(jat)

        dsi0d =   dpsi0(xd)*ddi_sav(iat)
        dsi1d =   dpsi1(xd)
        dsi2d =   dpsi2(xd)*dd_sav(iat)

        dsi0md = -dpsi0(mxd)*ddi_sav(iat)
        dsi1md =  dpsi1(mxd)
        dsi2md = -dpsi2(mxd)*dd_sav(iat)

        ! second derivatives of the weight functions
        ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
        ddsi1t =   ddpsi1(xt)*dti_sav(jat)
        ddsi2t =   ddpsi2(xt)

        ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
        ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
        ddsi2mt =  ddpsi2(mxt)

        ! free energy
        free  = h5(iat,jat, &
                   si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                   si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
        ! derivative with respect to density
        df_d  = h5(iat,jat, &
                   si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                   dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)
        ! derivative with respect to temperature
        df_t  = h5(iat,jat, &
                   dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                   si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
        ! derivative with respect to temperature**2
        df_tt = h5(iat,jat, &
                   ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                   si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
        ! derivative with respect to temperature and density
        df_dt = h5(iat,jat, &
                   dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                   dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

        !.....Compute pressure and derivatives...................

        ! get the interpolation weight functions
        si0t   =  xpsi0(xt)
        si1t   =  xpsi1(xt)*dt_sav(jat)

        si0mt  =  xpsi0(mxt)
        si1mt  =  -xpsi1(mxt)*dt_sav(jat)

        si0d   =  xpsi0(xd)
        si1d   =  xpsi1(xd)*dd_sav(iat)

        si0md  =  xpsi0(mxd)
        si1md  =  -xpsi1(mxd)*dd_sav(iat)

        ! look in the pressure derivative only once
        fi(1)  = dpdf(iat,jat)
        fi(2)  = dpdf(iat+1,jat)
        fi(3)  = dpdf(iat,jat+1)
        fi(4)  = dpdf(iat+1,jat+1)
        fi(5)  = dpdft(iat,jat)
        fi(6)  = dpdft(iat+1,jat)
        fi(7)  = dpdft(iat,jat+1)
        fi(8)  = dpdft(iat+1,jat+1)
        fi(9)  = dpdfd(iat,jat)
        fi(10) = dpdfd(iat+1,jat)
        fi(11) = dpdfd(iat,jat+1)
        fi(12) = dpdfd(iat+1,jat+1)
        fi(13) = dpdfdt(iat,jat)
        fi(14) = dpdfdt(iat+1,jat)
        fi(15) = dpdfdt(iat,jat+1)
        fi(16) = dpdfdt(iat+1,jat+1)

        ! pressure derivative with density
        dpepdd = h3(iat,jat, &
                    si0t,   si1t,   si0mt,   si1mt, &
                    si0d,   si1d,   si0md,   si1md)
        dpepdd = max(ye * dpepdd,1.0d-30)

        !..................................................................
        !.....Compute the desired e+-e- thermodynamic quantities...........

        ! dpepdd at high temperatures and low densities is below the
        ! floating point limit of the subtraction of two large terms.
        ! since dpresdd doesn't enter the maxwell relations at all, use the
        ! bicubic interpolation done above instead of the formally correct expression
        x       = din * din
        pele    = x * df_d
        dpepdt  = x * df_dt
!        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)

        x       = ye * ye
        sele    = -df_t * ye
        dsepdt  = -df_tt * ye
        dsepdd  = -df_dt * x

        eele    = ye*free + temp * sele
        deepdt  = temp * dsepdt
        deepdd  = x * df_d + temp * dsepdd

        !......................................................................
        !.....Coulomb section: account for charged interactions................
        !......................................................................

        ! uniform background corrections only from yakovlev & shalybkov 1989
        ! lami is the average ion seperation
        ! plasg is the plasma coupling parameter

        z        = forth * pi
        s        = z * xni
        dsdd     = z * dxnidd

        lami     = 1.0d0/s**third
        inv_lami = 1.0d0/lami
        z        = -third * lami
        lamidd   = z * dsdd/s

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        z        = -plasg * inv_lami
        plasgdd  = z * lamidd
        plasgdt  = -plasg*ktinv * kboltz

        ! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg .ge. 1.0) then
            x        = plasg**(0.25d0)
            y        = avo_real * ytot1 * kboltz
            ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
            pcoul    = third * den * ecoul
            scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
                       + d1 * (log(plasg) - 1.0d0) - e1)

            y        = avo_real*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
            decouldd = y * plasgdd
            decouldt = y * plasgdt + ecoul/temp

            y        = third * den
            dpcouldd = third * ecoul + y*decouldd
            dpcouldt = y * decouldt

            y        = -avo_real*kboltz/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
            dscouldd = y * plasgdd
            dscouldt = y * plasgdt

        ! yakovlev & shalybkov 1989 equations 102, 103, 104
        else if (plasg .lt. 1.0) then
            x        = plasg*sqrt(plasg)
            y        = plasg**b2
            z        = c2 * x - third * a2 * y
            pcoul    = -pion * z
            ecoul    = 3.0d0 * pcoul/den
            scoul    = -avo_real/abar*kboltz*(c2*x -a2*(b2-1.0d0)/b2*y)

            s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
            dpcouldd = -dpiondd*z - pion*s*plasgdd
            dpcouldt = -dpiondt*z - pion*s*plasgdt

            s        = 3.0d0/den
            decouldd = s * dpcouldd - ecoul/den
            decouldt = s * dpcouldt

            s        = -avo_real*kboltz/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
            dscouldd = s * plasgdd
            dscouldt = s * plasgdt
        end if

        ! if Coulomb corrections cause a negative pressure, set them to zero
        x = prad + pion + pele + pcoul
        y = erad + eion + eele + ecoul
        z = srad + sion + sele + scoul
        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
            pcoul    = 0.0d0
            dpcouldd = 0.0d0
            dpcouldt = 0.0d0
            ecoul    = 0.0d0
            decouldd = 0.0d0
            decouldt = 0.0d0
            scoul    = 0.0d0
            dscouldd = 0.0d0
            dscouldt = 0.0d0
        end if

        !......................................................................
        !.....Sum all the contributions........................................
        !......................................................................

        ! sum all the gas components
        pgas    = pion + pele + pcoul
        egas    = eion + eele + ecoul
        sgas    = sion + sele + scoul

        dpgasdd = dpiondd + dpepdd + dpcouldd
        dpgasdt = dpiondt + dpepdt + dpcouldt

        degasdd = deiondd + deepdd + decouldd
        degasdt = deiondt + deepdt + decouldt

        dsgasdd = dsiondd + dsepdd + dscouldd
        dsgasdt = dsiondt + dsepdt + dscouldt

        ! add in radiation to get the total
        pres    = prad + pgas
        ener    = erad + egas
        entr    = srad + sgas

        dpresdd = dpraddd + dpgasdd
        dpresdt = dpraddt + dpgasdt

        denerdd = deraddd + degasdd
        denerdt = deraddt + degasdt

        dentrdd = dsraddd + dsgasdd
        dentrdt = dsraddt + dsgasdt

        !......................................................................
        !.....Compute specific heats and thermodynamic exponents...............
        !......................................................................

        ! compute the sound velocity
        zz    = pres*deni
        zzi   = den/pres
        chit  = temp/pres * dpresdt
        chid  = dpresdd*zzi
        cv    = denerdt
        x     = zz * chit/(temp * cv)
        gam1  = chit*x + chid
        z     = 1.0d0 + (ener + light2)*zzi
        sound = clite * sqrt(gam1/z)

        !......................................................................
        !.....Final checks and storing.........................................
        !......................................................................

        ! maxwell relations; each is zero if the consistency is perfect
        x   = den * den
        dse = temp*dentrdt/denerdt - 1.0d0
        dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
        dsp = -dentrdd*x/dpresdt - 1.0d0

        ! maxwel test
        !    dse_row(j)    = dse
        !    dpe_row(j)    = dpe
        !    dsp_row(j)    = dsp

        ! store outputs (new part, used in kNEC)
        px(j)    = pres
        ex(j)    = ener
        cs2x(j)  = sound * sound
        dpdtx(j) = dpresdt
        dedtx(j) = denerdt
        sx(j)    = entr
        pradx(j) = prad

      end do

    end subroutine helmholtz_eos

end module helmholtz
