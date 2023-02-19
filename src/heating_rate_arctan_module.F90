!***********************************************************************
!
!     module: heating_rate_module
!
!     This module contains a set of subroutines to read in fitted
!     heating rate, to interpolate them on a regular grid, and to
!     compute heating rate at selected times
!      
!     The heating rate is assumed to be in the form:
!              
!        eps(t) = B(1) * (1/2 - 1\pi*arctan((t-arctan_t0)/sigma))**B(2) for t <= arctan_t2
!        eps(t) = A(1) * t**(-A(2))   for t > arctan_t2
!
!     where:
!      
!        t    ...  time [s]
!        eps  ...  specific heating rate [erg/g/s]
!        B(2) ...  heating rate prefactor [erg/g/s]
!        A(1) ...  heating rate prefactor [erg/g/s]
!        B(2) ...  heating rate slope [-]
!        A(2) ...  heating rate slope [-]
!
!     A(1) and A(2) were interpolated through a best fit procedure on a
!     large sample of SkyNet trajectories. A(1) and A(2) are provided on
!     a regular 3D grid: A = A(tau,s,ye), where
!
!        tau ... expansion timescale [ms]
!        s   ... specific entropy [k_B/baryon]
!        ye  ... electron fraction [-]
!      
!     For the 3D grid variable,
!      
!        tau: 15 points between 1.36ms and 100ms
!        s  : 21 points between 1.82 k_B/baryon and 100 k_B/baryon
!        ye : 24 points between 0.02 and 0.48
!
!     The points in the tau and s grids are equally spaced in log10. The
!     ye points are equally spaced in linear scale. The tau, s and ye
!     grids are initialized by providing the minimal values and the
!     steps.
!      
!     A trilinear interpolation is performed for the A coefficients.
!
!     The arctan expression is taken from Korobkin's fit and it is done
!     in such a way that:
!      
!        * for t small, eps produces an almost flat behavior around B(1)
!        * eps has a small time behavior similar to what observed in
!          network trajectories, regulated by arctan_t0 and arctan_sigma0 (see below):
!          in particular, after the plateau, it has a drop and after that 
!          it entres a first power-law phase
!        * for t=arctan_t2, it connects almost continuosly to the (second,t>arctan_t2) 
!          power-law phase.
!          
!     The small time behavior is further regulated by a few parameters:
!       
!        - arctan_eps0     the asymptotic flat value, such that B(1) = arctan_eps0
!        - arctan_t0       determine the extension of the plateau phase
!        - arctan_sigma0   determines how much eps drops after the plateau,
!                   before entering the power-law phase
!       
!     If arctan_eps0 and arctan_sigma0 are assigned, B(2) is calculated to provide
!     a quasi-continuos expression for eps (it is assumed that arctan_t0 does
!     not matter for that)
!
!
!***********************************************************************

      module heating_rate_arctan_module
      use parameters, only: arctan_t2,arctan_eps0,arctan_sigma0,arctan_t0
      use physical_constants
      implicit none

      integer, parameter :: n_tau = 15               ! number of tau points
      integer, parameter :: n_s   = 21               ! number of entropy points
      integer, parameter :: n_ye  = 24               ! number of ye points

      real(8), dimension(n_tau) :: ltau_grid         ! log10(tau) grid
      real(8), dimension(n_tau) :: tau_grid          ! tau grid
      real(8), dimension(n_s)   :: ls_grid           ! log10(entropy) grid
      real(8), dimension(n_s)   :: s_grid            ! entropy grid
      real(8), dimension(n_ye)  :: ye_grid           ! ye grid

      real(8), parameter :: ltau_min = 1.2d0/9.d0    ! log10(tau_min)=0.133
      real(8), parameter :: ls_min   = 6.d0/2.3d1    ! log10(s_min)=0.2609
      real(8), parameter :: ye_min   = 2.d-2         ! ye_min=0.02

      real(8), parameter :: ltau_delta = 1.2d0/9.d0  ! delta(log10(tau))=0.133
      real(8), parameter :: ls_delta   = 2.d0/2.3d1  ! delta(log10(s))=0.087
      real(8), parameter :: ye_delta   = 2.d-2       ! delta(ye)=0.02

      real(8), parameter :: toll       = 1.d-4       ! grid tollerance
      real(8), parameter :: day2sec    = 8.64d4      ! day-to-sec conversion factor
      real(8), parameter :: sec2day    = 1.d0/day2sec ! sec-to-day conversion factor


!      real(8), parameter :: pi         = 2.d0*dacos(0.d0)  ! pi [-]
!      real(8), parameter :: oneoverpi  = 1.d0/pi     ! 1/pi [-]
      real(8), dimension(2,n_tau,n_s,n_ye) :: A      ! fitted coefficients

!.....early time behavior...............................................
!      real(8), parameter :: arctan_t2   = 0.1                 ! 0.1 day
!
!      ! used defined parameters !
!      real(8), parameter :: arctan_eps0       = 1.0d19        ! [erg/g/s]
!      real(8), parameter :: arctan_sigma0     = 0.11          ! heating fit parameter [s], Korobkin: 0.11
!      real(8), parameter :: arctan_t0         = 2.0           ! heating fit parameter [s], Korobkin: 1.3
!
      contains

!=======================================================================
!
!     subroutine: read_heating_table
!              
!     This subroutine reads in the table with the arrays of the
!     independent variables and with the interpolation coefficients
!
!=======================================================================

      subroutine read_heating_table_arctan(filename)

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

      real(8), dimension(2) :: A_tmp

      open(unit=11,file=trim(filename),status='old')

!.....skip the file header..............................................
      do i=1,7
        read(11,*)
      end do

!.....build the tau, s and ye grids.....................................
      do i=1,n_tau
        ltau_grid(i) = ltau_min + ltau_delta*dble(i-1)   ! log10 scale
      end do
      do i=1,n_s
        ls_grid(i) = ls_min + ls_delta*dble(i-1)         ! log10 scale
      end do
      do i=1,n_ye
        ye_grid(i) = ye_min + ye_delta*dble(i-1)         ! linear scale
      end do

!.....compute the tau and entropy grids also in linear scale      
      tau_grid = 10.d0**ltau_grid
      s_grid = 10.d0**ls_grid

!.....read in the full table............................................
      do i=1,n_tau
        do j=1,n_s
          do k=1,n_ye

            read(11,12)tau_tmp,s_tmp,ye_tmp,A_tmp
            A(:,i,j,k) = A_tmp

!...........check if the tau grid is compatible with the table tau......
            if (abs((tau_tmp-tau_grid(i))/tau_tmp).gt.toll) then
              write(6,*)'Compatibility problem for tau!'
              write(6,*)'table tau',tau_tmp
              write(6,*)'grid tau',tau_grid(i)
              write(6,*)'grid position',i
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
            if (abs((ye_tmp-ye_grid(k))/ye_tmp).gt.toll) then
              write(6,*)'Compatibility problem for ye!'
              write(6,*)'table ye',ye_tmp
              write(6,*)'grid ye',ye_grid(k)
              stop  
            end if

          end do
        end do
      end do
12    format(5es12.4)

      close(11)

      return

      end subroutine read_heating_table_arctan

!=======================================================================

!=======================================================================
!
!     subroutine: set_tau_index_arctan
!              
!     This subroutine computes the index of the grid entry immediately
!     below the expansion timescale tau. It additionally provides the   
!     normalized distance from that grid point.
!
!=======================================================================

      subroutine set_tau_index_arctan(tau,itau,d_tau)

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

!.....take the log10 of the expansion timescale
      tmp = log10(tau)

      if (tmp.le.ltau_grid(1)) then
        itau  = 1
        d_tau = 0.d0
      else if (tmp.ge.ltau_grid(n_tau)) then
        itau  = n_tau-1
        d_tau = 1.d0
      else
        itau = 1 + floor((tmp-ltau_min)/ltau_delta)
        d_tau = (tmp-ltau_grid(itau))/ltau_delta
      end if

      return

      end subroutine set_tau_index_arctan

!=======================================================================

!=======================================================================
!
!     subroutine: set_s_index_arctan
!              
!     This subroutine computes the index of the grid entry immediately
!     below the entropy s. It additionally provides the normalized 
!     distance from that grid point.
!
!=======================================================================

      subroutine set_s_index_arctan(s,is,d_s)

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

!.....take the log10 of the entropy
      tmp = log10(s)

      if (tmp.le.ls_grid(1)) then
        is  = 1
        d_s = 0.d0
      else if (tmp.ge.ls_grid(n_s)) then
        is  = n_s-1
        d_s = 1.d0
      else
        is = 1 + floor((tmp-ls_min)/ls_delta)
        d_s = (tmp-ls_grid(is))/ls_delta
      end if

      return

      end subroutine set_s_index_arctan

!=======================================================================

!=======================================================================
!
!     subroutine: set_ye_index_arctan
!              
!     This subroutine computes the index of the grid entry immediately
!     below the electron fraction ye. It additionally provides the   
!     normalized distance from that grid point.
!
!=======================================================================

      subroutine set_ye_index_arctan(ye,iye,d_ye)

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
        iye = 1 + floor((tmp-ye_min)/ye_delta)
        d_ye = (tmp-ye_grid(iye))/ye_delta
      end if

      return

      end subroutine set_ye_index_arctan

!=======================================================================

!=======================================================================
!
!     subroutine: interp_heating_coeff
!              
!     This subroutine interpolates the heating rate coefficients for
!     a set of specific conditions (tau,s,ye)
!
!=======================================================================

      subroutine interp_heating_coeff_arctan(tau,s,ye,A_int,B_int)

      implicit none

      real(8), intent(in)  :: tau
      real(8), intent(in)  :: s
      real(8), intent(in)  :: ye
      real(8), dimension(2), intent(out) :: A_int
      real(8), dimension(2), intent(out) :: B_int

!-----------------------------------------------------------------------
!
!     input:
!     tau ... expansion timescale [ms]
!     s   ... specific entropy [kB/baryon]
!     ye  ... electron fraction [-]
!     t   ... time [s]
!
!     output:
!     A  ... heating rate coefficients [erg/s/g or -]
!     B  ... heating rate coefficients [erg/s/g or -]
!      
!-----------------------------------------------------------------------


      integer :: itau,is,iye
      real(8) :: d_s,d_ye,d_tau
      real(8) :: onemdtau,onemds,onemdye

!.....find the relevant indexes
      call set_tau_index_arctan(tau,itau,d_tau)
      call set_s_index_arctan(s,is,d_s)
      call set_ye_index_arctan(ye,iye,d_ye)

      !write(6,*)tau,tau_grid(itau),tau_grid(itau+1),d_tau
      !write(6,*)s,s_grid(is),s_grid(is+1),d_s
      !write(6,*)ye,ye_grid(iye),ye_grid(iye+1),d_ye

      onemdtau = 1.d0-d_tau
      onemds   = 1.d0-d_s
      onemdye  = 1.d0-d_ye

      A_int = onemdye*(onemds*(onemdtau * A(:,itau  ,is  ,iye  )        &
       +                          d_tau * A(:,itau+1,is  ,iye  ) )      &
       +                  d_s*(onemdtau * A(:,itau  ,is+1,iye  )        &
       +                          d_tau * A(:,itau+1,is+1,iye  ) ) )    &
       +         d_ye*(onemds*(onemdtau * A(:,itau  ,is  ,iye+1)        &
       +                          d_tau * A(:,itau+1,is  ,iye+1) )      &
       +                  d_s*(onemdtau * A(:,itau  ,is+1,iye+1)        &
       +                          d_tau * A(:,itau+1,is+1,iye+1) ) )

      B_int(1) = arctan_eps0
      B_int(2) = (A_int(2)*dlog(arctan_t2)-dlog(A_int(1)/arctan_eps0))/dlog(arctan_t2*day2sec*pi/arctan_sigma0)

      return

      end subroutine interp_heating_coeff_arctan

!=======================================================================

!=======================================================================
!
!     subroutine: calc_heating_rate
!              
!     This subroutine computes the heating rate at a specific time, for
!     a set of specific conditions (tau,s,ye)
!
!=======================================================================

      subroutine calc_heating_rate_arctan(tau,s,ye,t,eps)

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

      real(8), dimension(2) :: A_int
      real(8), dimension(2) :: B_int

      real(8) :: t_day

      call interp_heating_coeff_arctan(tau,s,ye,A_int,B_int)

!.....compute the heating rate..........................................
      t_day = t*sec2day
      if (t_day.ge.arctan_t2) then
          eps = A_int(1)*t_day**(-A_int(2))
      else
          eps = B_int(1)*(0.5-oneoverpi*datan((t-arctan_t0)/arctan_sigma0))**B_int(2)
      end if

      return

      end subroutine calc_heating_rate_arctan

!=======================================================================
!=======================================================================
!
!     subroutine: calc_heating_rate_t_array
!              
!     This subroutine computes the heating rate for a set of specific 
!     conditions (tau,s,ye) along an array of times of size ntime
!
!=======================================================================

      subroutine calc_heating_rate_t_array_arctan(ntime,tau,s,ye,t_ar,eps_ar)

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
!     ntime  ... number of times in time array [-]
!     tau    ... expansion timescale [ms]
!     s      ... specific entropy [kB/baryon]
!     ye     ... electron fraction [-]
!     t_ar   ... time array [s]
!
!     output:
!     eps_ar  ... heating rate array [erg/g/s]
!      
!-----------------------------------------------------------------------

      real(8) :: arctan_eps0
      real(8), dimension(2) :: A_int,B_int

      real(8), dimension(ntime) :: t_day_ar

      call interp_heating_coeff_arctan(tau,s,ye,A_int,B_int)

!.....compute the heating rate..........................................
      t_day_ar = t_ar*sec2day
      where (t_day_ar.ge.arctan_t2)
          eps_ar = A_int(1)*t_day_ar**(-A_int(2))
      elsewhere
          eps_ar = B_int(1)*(0.5-oneoverpi*datan((t_ar-arctan_t0)/arctan_sigma0))**B_int(2)
      end where

      return

      end subroutine calc_heating_rate_t_array_arctan

!=======================================================================

      end module heating_rate_arctan_module

!***********************************************************************
