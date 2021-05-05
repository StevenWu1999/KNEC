!***********************************************************************
!
!     module: heating_rate_LR15_module
!
!     This module contains a set of subroutines to read in fitted
!     heating rate provided by Lippuner and Roberts 2015 (publicly 
!     available data), to interpolate them on a regular grid, and to
!     compute heating rate at selected times
!      
!     The heating rate is assumed to be in the form:
!              
!              eps(t) = A(1) * t**(-A(2))   + 
!                       A(3) * exp(-t/A(4)) +
!                       A(5) * exp(-t/A(6)) +
!                       A(7) * exp(-t/A(8))
!
!     It is important to notice that the fitting formula is not always
!     the same: according to the grid point, up to three of the four
!     terms in the above formula are used, i.e. sometimes one term, 
!     sometimes two, sometimes three.      
!     The actual fitting formula is the one that provides the best fit.
!
!     In the previous expression:
!      
!        t    ...  time [days]
!        eps  ...  specific heating rate [erg/g/s]
!        A(1) ...  power-law prefactor [erg/g/s]
!        A(2) ...  power-law rate slope [-]
!        A(3) ...  first exponential prefactor [erg/g/s]
!        A(4) ...  first exponential argument factor [-]
!        A(5) ...  first exponential prefactor [erg/g/s]
!        A(6) ...  first exponential argument factor [-]
!        A(7) ...  first exponential prefactor [erg/g/s]
!        A(8) ...  first exponential argument factor [-]
!
!     A(:) were interpolated through a best fit procedure on a
!     large sample of SkyNet trajectories. They are provided on
!     a 3D grid: A = A(ye,s,tau), where
!
!        ye  ... electron fraction [-]
!        s   ... specific entropy [k_B/baryon]
!        tau ... expansion timescale [ms]
!      
!     For the 3D grid variable,
!      
!        ye : 17 points between 0.01 and 0.50
!        s  : 17 points between 1.0 k_B/baryon and 100 k_B/baryon
!        tau: 17 points between 0.1ms and 500ms
!
!     The points in the tau and s grids are approximatelly equally 
!     spaced in log10. The ye points are equally spaced (delta=0.03) 
!     in linear scale, with a jump between 0.25 and 0.29 
!      
!     The different subroutines read in the fitted coefficeint and
!     provide interpolated rates for a specific triplet of tau, s and 
!     ye, and at a given time (or an array of times).
!     A trilinear interpolation is done on the log10 of the heating rate
!     on the ye, s, tau space.      
!     
!
!***********************************************************************

      module heating_rate_LR15_module

      implicit none

      integer, parameter :: n_ye  = 17               ! number of ye points
      integer, parameter :: n_s   = 17               ! number of entropy points
      integer, parameter :: n_tau = 17               ! number of tau points

      real(8), dimension(n_tau) :: ltau_grid         ! log10(tau) grid
      real(8), dimension(n_tau) :: tau_grid          ! tau grid
      real(8), dimension(n_s)   :: ls_grid           ! log10(entropy) grid
      real(8), dimension(n_s)   :: s_grid            ! entropy grid
      real(8), dimension(n_ye)  :: ye_grid           ! ye grid

      real(8), parameter :: ltau_min   = -1.d0       ! min log10(tau)
      real(8), parameter :: ltau_delta = 0.231185625 ! delta log10(tau)

      real(8), parameter :: ls_min   = 0.d0          ! min log10(s)
      real(8), parameter :: ls_delta = 0.125         ! delta log10(s)

      real(8), parameter :: ye_min   = 0.01          ! min log10(s)
      real(8), parameter :: ye_delta = 0.03          ! delta log10(s)

      real(8), parameter :: toll       = 1.d-4       ! grid tollerance
      real(8), parameter :: sec2day    = 1.d0/8.64e4 ! sec-to-day conversion factor

      real(8), dimension(8,n_tau,n_s,n_ye) :: A      ! fitted coefficients

      contains

!=======================================================================
!
!     subroutine: read_heating_table
!              
!     This subroutine reads in the table with the arrays of the
!     independent variables and with the interpolation coefficients
!     Note that for LR15 data the variable order is (ye,s,tau).
!     However, they are stored as A(tau,s,ye)
!              
!=======================================================================

      subroutine read_heating_table_LR15(filename)

      implicit none

      character(80), intent(in) :: filename

!-----------------------------------------------------------------------
!
!     input:
!     filename ... path and name of the fitted table      
!
!-----------------------------------------------------------------------

      integer :: i,j,k
      real(8) :: tau_tmp,s_tmp,ye_tmp

      real(8), dimension(8) :: A_tmp

      open(unit=11,file=trim(filename),status='old')

!.....skip the file header..............................................
      do i=1,13
        read(11,*)
      end do

!.....build the tau, s and ye grids.....................................
      tau_grid=(/ 0.1, 0.17, 0.29, 0.49, 0.84, 1.4, 2.4, 4.2, 7.1, 12., &
                  21., 35. , 59., 100., 170., 290., 500. /)
      ltau_grid = dlog10(tau_grid)


      s_grid=(/ 1.0, 1.3, 1.8, 2.4, 3.2, 4.2, 5.6, 7.5, 10., 13., 18.,  &
              24., 32., 42., 56., 75., 100. /)
      ls_grid = dlog10(s_grid)

      do i=1,n_ye
        ye_grid(i) = ye_min + ye_delta*dble(i-1)         ! linear scale
      end do
!.....correct the ye grid...............................................
      ye_grid(10:) = ye_grid(10:)+1.d-2

!.....read in the full table............................................
      do i=1,n_ye
        do j=1,n_s
          do k=1,n_tau

            read(11,*)ye_tmp,s_tmp,tau_tmp,A_tmp
            A(:,k,j,i) = A_tmp

!...........check if the tau grid is compatible with the table tau......
            if (abs((tau_tmp-tau_grid(k))/tau_tmp).gt.toll) then
              write(6,*)'Compatibility problem for tau!'
              write(6,*)'table tau',tau_tmp
              write(6,*)'grid tau',tau_grid(k)
              stop  
            end if

!...........check if the s grid is compatible with the table s..........
            if (abs((s_tmp-s_grid(j))/s_tmp).gt.toll) then
              write(6,*)'Compatibility problem for entropy!'
              write(6,*)'table s',s_tmp
              write(6,*)'grid s',s_grid(j)
              stop  
            end if

!...........check of the ye grid is compatible with the table ye........
            if (abs((ye_tmp-ye_grid(i))/ye_tmp).gt.toll) then
              write(6,*)'Compatibility problem for ye!'
              write(6,*)'table ye',ye_tmp
              write(6,*)'grid ye',ye_grid(i)
              stop  
            end if

          end do
        end do
      end do
12    format(5es12.4)

      close(11)

      return

      end subroutine read_heating_table_LR15

!=======================================================================

!=======================================================================
!
!     subroutine: set_tau_index
!              
!     This subroutine computes the index of the grid entry immediately
!     below the expansion timescale tau. It additionally provides the   
!     normalized distance from that grid point.
!
!=======================================================================

      subroutine set_tau_index_LR15(tau,itau,d_tau)

      implicit none

      real(8), intent(in)  :: tau
      integer, intent(out) :: itau
      real(8), intent(out) :: d_tau

!-----------------------------------------------------------------------
!
!     input:
!     tau ... expansion timescale [ms]
!
!     output:
!     itau  ... index of the closest smaller grid entry on the tau array [-]
!     d_tau ... normalized distance of tau from the itau point [-]
!      
!-----------------------------------------------------------------------

      real(8) :: tmp

!.....take the log of the expansion timescale...........................
      tmp = log10(tau)

      if (tmp.le.ltau_grid(1)) then
        itau  = 1
        d_tau = 0.d0
      else if (tmp.ge.ltau_grid(n_tau)) then
        itau  = n_tau-1
        d_tau = 1.d0
      else
!.......guess the time index............................................
        itau = 1 + floor((tmp-ltau_min)/ltau_delta)

        if (ltau_grid(itau).lt.tmp) then
          continue
        else
          itau = itau -1 
          if (ltau_grid(itau).lt.tmp) then
             continue
          else
             write(6,*)'Something wrong with tau grid!'
             stop
          end if
        end if
        d_tau = (tmp-ltau_grid(itau))/(ltau_grid(itau+1)-ltau_grid(itau))
      end if

      return

      end subroutine set_tau_index_LR15

!=======================================================================

!=======================================================================
!
!     subroutine: set_s_index
!              
!     This subroutine computes the index of the grid entry immediately
!     below the entropy s. It additionally provides the normalized 
!     distance from that grid point.
!
!=======================================================================

      subroutine set_s_index_LR15(s,is,d_s)

      implicit none

      real(8), intent(in)  :: s
      integer, intent(out) :: is
      real(8), intent(out) :: d_s

!-----------------------------------------------------------------------
!
!     input:
!     s ... entropy [kb/baryon]
!
!     output:
!     is  ... index of the closest smaller grid entry on the s array [-]
!     d_s ... normalized distance of s from the is grid points [-]
!      
!-----------------------------------------------------------------------

      real(8) :: tmp

!.....take the log10 of the entropy.....................................
      tmp = log10(s)

      if (tmp.le.ls_grid(1)) then
        is  = 1
        d_s = 0.d0
      else if (tmp.ge.ls_grid(n_s)) then
        is  = n_s-1
        d_s = 1.d0
      else
        is = 1 + floor((tmp-ls_min)/ls_delta)
        if (ls_grid(is).lt.tmp) then
          continue
        else
          is = is -1 
          if (ls_grid(is).lt.tmp) then
             continue
          else
             write(6,*)'Something wrong with s grid!'
             stop
          end if        
        end if        
        d_s = (tmp-ls_grid(is))/(ls_grid(is+1)-ls_grid(is))
      end if

      return

      end subroutine set_s_index_LR15

!=======================================================================

!=======================================================================
!
!     subroutine: set_ye_index_LR15
!              
!     This subroutine computes the index of the grid entry immediately
!     below the electron fraction ye. It additionally provides the   
!     normalized distance from that grid point.
!
!=======================================================================

      subroutine set_ye_index_LR15(ye,iye,d_ye)

      implicit none

      real(8), intent(in)  :: ye
      integer, intent(out) :: iye
      real(8), intent(out) :: d_ye

!-----------------------------------------------------------------------
!
!     input:
!     ye ... electron fraction [-]
!
!     output:
!     iye  ... index of the closest smaller grid entry on the ye array [-]
!     d_ye ... normalized distance of ye from the iye grid points [-]
!      
!-----------------------------------------------------------------------

      real(8) :: tmp

      tmp = ye
      if (tmp.le.ye_grid(1)) then
        iye  = 1
        d_ye = 0.d0
      else if (tmp.ge.ye_grid(n_ye)) then
        iye  = n_ye-1
        d_ye = 1.d0
      else
        if (tmp.lt.0.25) then
          iye = 1 + floor((tmp-ye_min)/ye_delta)
          d_ye = (tmp-ye_grid(iye))/(ye_grid(iye+1)-ye_grid(iye))
        else if (tmp.gt.0.29) then 
          iye = 10 + floor((tmp-0.29)/ye_delta)
          d_ye = (tmp-ye_grid(iye))/(ye_grid(iye+1)-ye_grid(iye))
        else
          iye = 9
          d_ye = (tmp-ye_grid(iye))/0.04
        end if
      end if

      return

      end subroutine set_ye_index_LR15

!=======================================================================

!=======================================================================
!
!     subroutine: calc_eps_cube_LR15
!              
!     This subroutine computes the log10 of the heating rates on a cube 
!     in the LR15 grid specified by the indexes of the vertex (1,1,1),
!     where the general vertix is given by (i,j,k) with i,j,k=1,2 
!
!=======================================================================

      subroutine calc_eps_cube_LR15(itau,is,iye,t_day,leps)

      implicit none

      integer, intent(in) :: itau
      integer, intent(in) :: is
      integer, intent(in) :: iye
      real(8), intent(in) :: t_day
      real(8), dimension(2,2,2), intent(out) :: leps

!-----------------------------------------------------------------------
!
!     input:
!     itau  ... expansion timescale index [-]
!     is    ... specific entropy index    [-]
!     iye   ... electron fraction index   [-]
!     t_day ... time                      [days]
!
!     output:
!     leps  ... log10 of the heating rates in the first neighbours [erg/s/g]
!      
!-----------------------------------------------------------------------

      integer :: i,j,k
      integer :: ii,jj,kk

      do k=1,2
        kk = iye-1+k
        do j=1,2
          jj = is-1+j
          do i=1,2
            ii = itau-1+i
            leps(i,j,k) = A(1,ii,jj,kk)*t_day**(-A(2,ii,jj,kk)) +       &
                          A(3,ii,jj,kk)*dexp(-t_day/A(4,ii,jj,kk)) +    &
                          A(5,ii,jj,kk)*dexp(-t_day/A(6,ii,jj,kk)) +    &
                          A(7,ii,jj,kk)*dexp(-t_day/A(8,ii,jj,kk))

          end do
        end do
      end do

!.....take the log 10 of the heating rates..............................
      leps = dlog10(leps)

      return

      end subroutine calc_eps_cube_LR15

!=======================================================================

!=======================================================================
!
!     subroutine: calc_heating_rate_LR15
!              
!     This subroutine computes the heating rate at a specific time, for
!     a set of specific conditions (tau,s,ye)
!
!=======================================================================

      subroutine calc_heating_rate_LR15(tau,s,ye,t,eps)

      implicit none

      real(8), intent(in)  :: tau
      real(8), intent(in)  :: s
      real(8), intent(in)  :: ye
      real(8), intent(in)  :: t
      real(8), intent(out) :: eps

!-----------------------------------------------------------------------
!
!     input:
!     tau ... expansion timescale [ms]
!     s   ... specific entropy [kB/baryon]
!     ye  ... electron fraction [-]
!     t   ... time [s]
!
!     output:
!     eps  ... heating rate [erg/g/s]
!      
!-----------------------------------------------------------------------


      integer :: itau,is,iye
      real(8) :: d_tau,d_s,d_ye
      real(8) :: onemdtau,onemds,onemdye

      real(8) :: t_day,eps_int
      real(8), dimension(2,2,2) :: leps


!.....find the relevant indexes
      call set_tau_index_LR15(tau,itau,d_tau)
      call set_s_index_LR15(s,is,d_s)
      call set_ye_index_LR15(ye,iye,d_ye)

      onemdtau = 1.d0 - d_tau
      onemds   = 1.d0 - d_s
      onemdye  = 1.d0 - d_ye

!      write(6,*)'itau',itau
!      write(6,*)'is',is
!      write(6,*)'iye',iye
!
!      write(6,*)
!      write(6,*)'dtau',d_tau
!      write(6,*)'ds',d_s
!      write(6,*)'dye',d_ye
!
!      write(6,*)
!      write(6,*)'1-dtau',onemdtau
!      write(6,*)'1-ds',onemds
!      write(6,*)'1-dye',onemdye

!.....convert the time in days..........................................
      t_day = t*sec2day
     
!.....compute the heating rates on the first neighbours.................
      call calc_eps_cube_LR15(itau,is,iye,t_day,leps)
      


!.....interpolate linearly the log of the heating rates.................
      eps_int = onemdtau*(onemds*(onemdye * leps(1,1,1)        &
       +                             d_ye * leps(1,1,2) )      &
       +                     d_s*(onemdye * leps(1,2,1)        &
       +                             d_ye * leps(1,2,2) ) )    &
       +           d_tau*(onemds*(onemdye * leps(2,1,1)        &
       +                             d_ye * leps(2,1,2) )      &
       +                     d_s*(onemdye * leps(2,2,1)        &
       +                             d_ye * leps(2,2,2) ) )

      eps = 10.d0**(eps_int)

      return

      end subroutine calc_heating_rate_LR15

!=======================================================================

!=======================================================================
!
!     subroutine: calc_heating_rate_t_array_LR15
!              
!     This subroutine computes the heating rate for a set of specific 
!     conditions (tau,s,ye) along an array of times of size ntime
!
!=======================================================================

      subroutine calc_heating_rate_t_array_LR15(ntime,tau,s,ye,t_ar,eps_ar)

      implicit none

      integer, intent(in)                    :: ntime
      real(8), intent(in)                    :: tau
      real(8), intent(in)                    :: s
      real(8), intent(in)                    :: ye
      real(8), dimension(ntime), intent(in)  :: t_ar
      real(8), dimension(ntime), intent(out) :: eps_ar

!-----------------------------------------------------------------------
!
!     input:
!     ntime  ... number of time in time array [-]
!     tau    ... expansion timescale [ms]
!     s      ... specific entropy [kB/baryon]
!     ye     ... electron fraction [-]
!     t_ar   ... time array [s]
!
!     output:
!     eps_ar  ... heating rate array [erg/g/s]
!      
!-----------------------------------------------------------------------

      real(8), dimension(ntime) :: t_day_ar
      real(8), dimension(2,2,2) :: leps
      integer :: i
      integer :: itau,is,iye
      real(8) :: d_tau,d_s,d_ye
      real(8) :: onemdtau,onemds,onemdye

!.....compute the time in days..........................................
      t_day_ar = t_ar*sec2day

!.....find the relevant indexes.........................................
      call set_tau_index_LR15(tau,itau,d_tau)
      call set_s_index_LR15(s,is,d_s)
      call set_ye_index_LR15(ye,iye,d_ye)

      onemdtau = 1.d0 - d_tau
      onemds   = 1.d0 - d_s
      onemdye  = 1.d0 - d_ye

      do i=1,ntime

!.......compute the heating rates on the first neighbours..............
        call calc_eps_cube_LR15(itau,is,iye,t_day_ar(i),leps)

!.......interpolate linearly the log of the heating rates...............
        eps_ar(i) = onemdtau*(onemds*(onemdye * leps(1,1,1)             &
       +                                 d_ye * leps(1,1,2) )           &
       +                         d_s*(onemdye * leps(1,2,1)             &
       +                                 d_ye * leps(1,2,2) ) )         &
       +               d_tau*(onemds*(onemdye * leps(2,1,1)             &
       +                                 d_ye * leps(2,1,2) )           &
       +                         d_s*(onemdye * leps(2,2,1)             &
       +                                 d_ye * leps(2,2,2) ) )

      end do

      eps_ar = 10.d0**(eps_ar)

      return

      end subroutine calc_heating_rate_t_array_LR15

!=======================================================================

      end module heating_rate_LR15_module

!***********************************************************************
