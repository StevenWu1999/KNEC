!***********************************************************************
!
!     module: heating_rate_module
!
!     This module contains a set of subroutines: 
!       to read in fitted heating rate, 
!       to interpolate them on a regular grid, 
!       to compute heating rate at selected times
!      
!     The heating rate is assumed to be in the form:
!              
!              eps(t) = B(3)*(0.5-ARCTAN((t-B(5))/B(6))/pi)**(B(4))  for t <= t1
!              eps(t) = B(1) * t**(-B(2))                            for t > t2
!              eps(t) = a continous interpolation in between         for t1 < t <= t2
!
!     To improve accuracy, some parameters are taken in log scale:
!      
!              eps(t) = 10**(A(3)+A(4)*log10(0.5-ARCTAN((t-A(5))/A(6))/pi))  for t <= t1
!                  A(3) = log10(B(3))
!                  A(4) = B(4)
!                  A(5) = B(5)
!                  A(6) = B(6)
!      
!              eps(t) = 10**(A(1) - A(2) * log10(t))                 for t > t2
!                  A(1) = log10(B(1))
!                  A(2) = B(2)
!
!     where:
!      
!        t    ...  time [s]
!        eps  ...  specific heating rate [erg/g/s]
!        A(1) ...  log10 of the heating rate prefactor at late times [erg/g/s]
!        A(2) ...  heating rate slope at late times [-]
!        A(3) ...  log10 of the heating rate prefactor at early times [erg/g/s]
!        A(4) ...  time that set the early time plateau [s]
!        A(5) ...  time that set the transition time plateau --> power-law [s]
!        A(6) ...  asymptotic heating rate slope at early times [-]
!
!     A(1:6) were interpolated through a best fit procedure on a large 
!     sample of SkyNet trajectories 
!
!     The fit was done assuming the time in days. Since we prefer here
!     to work with seconds, the interpolated coefficients must be
!     corrected:
!
!     A(1) --> log10( B(1) * day2sec**(B(2)) )
!     A(3) --> log10( B(3) * day2sec**(B(4)) )
!
!     A(1:6) are provided on a regular 3D grid: A = A(tau,s,ye), where
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
!     t1 and t2 have been fixed for the interpolation procedure, based
!     on preliminary deta inspection 
!     The smooth transition is done by:
!
!     x = (t-t2)/(t1-t2)
!     eps1 = eps for t < t1      
!     eps2 = eps for t > t2      
!     eps = 10**(log10(eps1)*(1-x) + log10(eps2)*x)
!      
!***********************************************************************

      module heating_Apr2_module
      use physical_constants
      implicit none

!      real(8), parameter :: pi         = 2.d0*dacos(0.d0)  ! pi   We already have pi in physical_constants
      real(8), parameter :: day2sec    = 8.64e4            ! day-to-sec conversion factor
      real(8), parameter :: lday2sec   = dlog10(day2sec)   ! log10 of day-to-sec conversion factor
      real(8), parameter :: sec2day    = 1.d0/8.64e4       ! sec-to-day conversion factor

      ! (tau,s,ye) grid parameters

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


      real(8), dimension(6,n_tau,n_s,n_ye) :: A      ! fitted coefficients

!.....transition times..................................................
      real(8), parameter :: t1   = 1.0d3  !1.d+3             ! end of early time [s]
      real(8), parameter :: t2   = 4.0d4  !4.d+4             ! beginning of late time [s]

      contains

!=======================================================================
!
!     subroutine: read_heating_table_Apr2
!              
!     This subroutine reads in the table with the arrays of the
!     independent variables and with the interpolation coefficients
!
!=======================================================================

      subroutine read_heating_table_Apr2(filename,filename_early)

      implicit none

      character(80), intent(in) :: filename
      character(80), intent(in) :: filename_early

!-----------------------------------------------------------------------
!
!     input:
!     filename ... path and name of the fitted table at late time 
!     filename_early ... path and name of the fitted table at early time
!
!-----------------------------------------------------------------------

      integer :: i,j,k,line
      real(8) :: tau_tmp,s_tmp,ye_tmp
      real(8) :: tau_tmp_early,s_tmp_early,ye_tmp_early

      real(8), dimension(2) :: A_tmp
      real(8), dimension(4) :: A_tmp_early

      open(unit=11,file=trim(filename),status='old')
      open(unit=21,file=trim(filename_early),status='old')

!.....skip the file header..............................................
      do i=1,7
        read(11,*)
      end do

      do i=1,9
        read(21,*)
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

      !line = 0
!.....read in the full table............................................
      do i=1,n_tau
        do j=1,n_s
          do k=1,n_ye

            !line = line + 1
            !write(6,*)line

            ! read in the coefficient for the late time behavior
            read(11,12)tau_tmp,s_tmp,ye_tmp,A_tmp
            ! correct for the time units
            A_tmp(1) = dlog10(A_tmp(1))+A_tmp(2)*lday2sec
            ! store the coefficients
            A(1:2,i,j,k) = A_tmp

            ! read in the coefficient for the early time behavior
            read(21,22)tau_tmp_early,s_tmp_early,ye_tmp_early,A_tmp_early
            ! correct for the time units
            A_tmp_early(1) = dlog10(A_tmp_early(1))+A_tmp_early(2)*lday2sec
            ! store the coefficients
            A(3:6,i,j,k) = A_tmp_early

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

!...........check if the two tables are compatible among them...........
            if (abs((ye_tmp-ye_tmp_early)/ye_tmp).gt.toll) then
              write(6,*)'Compatibility problem between the two tables'
              write(6,*)'early table ye',ye_tmp_early
              write(6,*)'late table ye',ye_tmp
              stop  
            end if

!...........check if the two tables are compatible among them...........
            if (abs((s_tmp-s_tmp_early)/s_tmp).gt.toll) then
              write(6,*)'Compatibility problem between the two tables'
              write(6,*)'early table s',s_tmp_early
              write(6,*)'late table s',s_tmp
              stop  
            end if

!...........check if the two tables are compatible among them...........
            if (abs((tau_tmp-tau_tmp_early)/tau_tmp).gt.toll) then
              write(6,*)'Compatibility problem between the two tables'
              write(6,*)'early table tau',tau_tmp_early
              write(6,*)'late table tau',tau_tmp
              stop  
            end if

          end do
        end do
      end do
12    format(5es12.4)
22    format(7es12.4)

      close(11)
      close(12)

      return

      end subroutine read_heating_table_Apr2

!=======================================================================

!=======================================================================
!
!     subroutine: set_tau_index_Apr2
!              
!     This subroutine computes the index of the grid entry immediately
!     below the expansion timescale tau. It additionally provides the   
!     normalized distance from that grid point.
!
!=======================================================================

      subroutine set_tau_index_Apr2(tau,itau,d_tau)

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

      end subroutine set_tau_index_Apr2

!=======================================================================

!=======================================================================
!
!     subroutine: set_tau_index_Apr2
!              
!     This subroutine computes the index of the grid entry immediately
!     below the entropy s. It additionally provides the normalized 
!     distance from that grid point.
!
!=======================================================================

      subroutine set_s_index_Apr2(s,is,d_s)

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

      end subroutine set_s_index_Apr2

!=======================================================================

!=======================================================================
!
!     subroutine: set_ye_index_Apr2
!              
!     This subroutine computes the index of the grid entry immediately
!     below the electron fraction ye. It additionally provides the   
!     normalized distance from that grid point.
!
!=======================================================================

      subroutine set_ye_index_Apr2(ye,iye,d_ye)

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

      end subroutine set_ye_index_Apr2

!=======================================================================

!=======================================================================
!
!     subroutine: interp_heating_coeff_Apr2
!              
!     This subroutine interpolates the heating rate coefficients for
!     a set of specific conditions (tau,s,ye)
!
!=======================================================================

      subroutine interp_heating_coeff_Apr2(tau,s,ye,A_int)

      implicit none

      real(8)              , intent(in)  :: tau
      real(8)              , intent(in)  :: s
      real(8)              , intent(in)  :: ye
      real(8), dimension(6), intent(out) :: A_int

!-----------------------------------------------------------------------
!
!     input:
!     tau ... expansion timescale [ms]
!     s   ... specific entropy [kB/baryon]
!     ye  ... electron fraction [-]
!     t   ... time [s]
!
!     output:
!     A  ... heating rate coefficients [erg/s/g or s or -]
!      
!-----------------------------------------------------------------------


      integer :: itau,is,iye
      real(8) :: d_s,d_ye,d_tau
      real(8) :: onemdtau,onemds,onemdye

!.....find the relevant indexes
      call set_tau_index_Apr2(tau,itau,d_tau)
      call set_s_index_Apr2(s,is,d_s)
      call set_ye_index_Apr2(ye,iye,d_ye)

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

      return

      end subroutine interp_heating_coeff_Apr2

!=======================================================================

!=======================================================================
!
!     subroutine: calc_heating_rate_Apr2
!              
!     This subroutine computes the heating rate at a specific time, for
!     a set of specific conditions (tau,s,ye)
!
!=======================================================================

      subroutine calc_heating_rate_Apr2(tau,s,ye,t,eps)

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

      real(8), dimension(6) :: A_int

      real(8) :: x
      real(8) :: eps1,eps2

      call interp_heating_coeff_Apr2(tau,s,ye,A_int)

!.....compute the heating rate..........................................
      if (t.ge.t2) then
          eps = A_int(1) - A_int(2)*dlog10(t)
      else if (t.le.t1) then
          eps = A_int(3)+A_int(4)*                                      &
                 dlog10(5.d-1-datan((t-A_int(5))/A_int(6))/pi)
      else
          eps1 = A_int(3)+A_int(4)*                                     &
                 dlog10(5.d-1-datan((t-A_int(5))/A_int(6))/pi)
          eps2 = A_int(1) - A_int(2)*dlog10(t)
          !x = (t-t1)/(t2-t1)
          x = (dlog10(t/t1))/(dlog10(t2/t1))  ! log(eps)-log(t) fits better than log(eps)-linear(t)
          eps = eps1*(1.d0-x) + eps2*x
      end if

      eps = 10.d0**(eps)

      return

      end subroutine calc_heating_rate_Apr2

!=======================================================================

!=======================================================================
!
!     subroutine: calc_heating_rate_Apr2_t_array
!              
!     This subroutine computes the heating rate for a set of specific 
!     conditions (tau,s,ye) along an array of times of size ntime
!
!=======================================================================

      subroutine calc_heating_rate_Apr2_t_array(ntime,tau,s,ye,t_ar,eps_ar)

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

      real(8), dimension(6)     :: A_int
      real(8), dimension(ntime) :: x
      real(8), dimension(ntime) :: eps1
      real(8), dimension(ntime) :: eps2

      call interp_heating_coeff_Apr2(tau,s,ye,A_int)

!.....compute the heating rate..........................................
      where (t_ar.ge.t2)
          eps_ar = A_int(1) - A_int(2)*dlog10(t_ar)
      else where (t_ar.le.t1)
          eps_ar = A_int(3)+A_int(4)*                                   &
                 dlog10(0.5d0-datan((t_ar-A_int(5))/A_int(6))/pi)
      elsewhere
          eps1 = A_int(3)+A_int(4)*                                     &
                 dlog10(5.d-1-datan((t_ar-A_int(5))/A_int(6))/pi)
          eps2 = A_int(1) - A_int(2)*dlog10(t_ar)
          !x = (t_ar-t1)/(t2-t1)
          x = (dlog10(t_ar/t1))/(dlog10(t2/t1))  ! log(eps)-log(t) fits better than log(eps)-linear(t)
          eps_ar = eps1*(1.d0-x) + eps2*x
      end where
      eps_ar = 10.d0**eps_ar


      return

      end subroutine calc_heating_rate_Apr2_t_array

!=======================================================================

      end module heating_Apr2_module

!***********************************************************************
