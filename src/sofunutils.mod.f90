module md_sofunutils
  !/////////////////////////////////////////////////////////////////////////
  ! Contains various utility functions
  !-------------------------------------------------------------------------
  use md_params_core

  implicit none

contains

  function dampen_variability( var, tau, var_memory ) result( out_memory )
    !/////////////////////////////////////////////////////////////////////////
    ! Calculates the updated variable accounting for a memory time scale tau.
    ! Following Eq. 5 in Makela et al. (2004) Tree Physiology 24, 369–376
    ! 
    ! d(var_memory) / dt = (1 / tau) * (var - var_memory)
    ! 
    !-------------------------------------------------------------------------
    ! arguments
    real, intent(in)    :: var           ! fast-varying variable
    real, intent(in)    :: tau           ! memory e-folding time scale (d)
    real, intent(in)    :: var_memory    ! damped (low-pass filtered) variable

    ! function return variable
    real :: out_memory

    ! local variable
    real :: dvar

    dvar = (1.0/tau) * (var - var_memory)
    out_memory = var_memory + dvar

  end function dampen_variability


  function running( presval, inow, window_length, method, prevval ) result( runningval )
    !/////////////////////////////////////////////////////////////////////////
    ! Returns running sum/average of a window of length 'window_length'
    ! and whose last element is the 'inow' item of 'presval'.
    ! If provided, 'prevval' contains the values for the previous year used for padding.
    ! Otherwise no padding is used.
    !-------------------------------------------------------------------------
    ! arguments
    real, dimension(ndayyear), intent(in) :: presval          ! Vector containing values for present year
    integer, intent(in) :: inow                               ! Index corresponding to "now" (day of year in present year). Negative values are allowed wihen prevval is provided.
    integer, intent(in) :: window_length                      ! Window length
    character(len=*), intent(in) :: method                    ! Either "sum" or "mean" for running sum or running mean
    real, dimension(ndayyear), optional, intent(in):: prevval ! Vector containing values for the previous year

    ! local variables
    real, dimension(2*ndayyear) :: valbuf
    integer :: idx_start, idx_end, effective_window_length

    ! function return variable
    real :: runningval

    ! Initialize indexes
    idx_end = ndayyear + inow
    idx_start = idx_end - window_length + 1

    ! Initialization of the buffer
    if (present(prevval)) then
      valbuf(1:ndayyear) = prevval
      effective_window_length = window_length
    else
      valbuf(1:ndayyear) = 0.0
      effective_window_length = MIN(window_length, inow)
    end if
    valbuf(ndayyear + 1 : 2 * ndayyear) = presval

    runningval = sum(valbuf(idx_start:idx_end))

    if (method == "mean") then
      if (effective_window_length <= 0) then
        ! Negative effective window lenght is not allowed. We set it 0 so that the computation of the mean will fail.
        effective_window_length = 0
      end if
      runningval = runningval / effective_window_length
    end if

  end function running

  subroutine aggregate(out, in, ratio)
    !////////////////////////////////////////////////////////////////
    ! Aggregate array 'in' with ratio 'ratio' using average scheme.
    !----------------------------------------------------------------
    ! arguments
    real, dimension(:), intent(inout) :: out      ! Output array (should have size of 'in' / ratio)
    real, dimension(:), intent(in) :: in          ! Input array
    integer, intent(in) :: ratio                   ! Sampling ratio

    ! local variables
    integer :: idx, idx_in_start, idx_in_end

    if (ratio == 1) then
      out(:) = in(:)
    else
      do idx = 1, SIZE(in)/ratio
        idx_in_start = (idx - 1) * ratio + 1
        idx_in_end = idx * ratio
        out(idx) = SUM(in(idx_in_start : idx_in_end))/ratio
      end do
    end if

  end subroutine aggregate

  function area( lat, dx, dy ) result( out_area )
    !////////////////////////////////////////////////////////////////
    ! Calculates grid cell area in m2 on a spherical Earth
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: lat      ! latitude (degrees N)
    real, intent(in) :: dx       ! resolution in longitude (degrees)
    real, intent(in) :: dy       ! resolution in latitude (degrees)

    ! function return variable
    real :: out_area             ! gridcell surface area (m2)

    ! local variables
    real, parameter :: r_earth = 6370000

    out_area = 4.0 * r_earth**2 * 0.5 * dx * pi / 180.0 * cos( abs(lat) * pi / 180.0 ) * sin( 0.5 * dy * pi / 180.0 )

  end function area


  function calc_patm( elv ) result( patm )
    !----------------------------------------------------------------   
    ! Calculates atmospheric pressure for a given elevation, assuming
    ! standard atmosphere at sea level (kPo)
    ! Ref:      Allen et al. (1998)
    ! This function is copied from SPLASH
    !----------------------------------------------------------------
    ! arguments
    real, intent(in) :: elv ! elevation above sea level, m

    ! function return value
    real ::  patm ! atmospheric pressure (Pa)

    patm = kPo * (1.0 - kL * elv / kTo) ** (kG * kMa * 1.e-3 / (kR * kL))

  end function calc_patm

  function calc_esat(T) result( out_esat ) ! pressure, Pa
    implicit none
    real :: out_esat
    real, intent(in) :: T ! degC

    out_esat=610.78*exp(17.27*T/(T+237.3))

  end function calc_esat

  function median( vec, len ) result( out )
    !--------------------------------------------------------------------
    ! This function receives an array vec of N entries, copies its value
    ! to a local array Temp(), sorts Temp() and computes the median.
    ! The returned value is of REAL type.
    ! Copied from https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/median.f90
    !--------------------------------------------------------------------
    real, dimension(:), intent(in)     :: vec
    integer, intent(in)                :: len

    ! function return variable
    real                               :: out

    ! local variables
    real, dimension(len)             :: tmp

    tmp(:) = vec(:)

    call  sort(tmp, len)                ! sort the copy

    if (mod(len,2) == 0) then           ! compute the median
      out = (tmp(len/2) + tmp(len/2+1)) / 2.0
    else
      out = tmp(len/2+1)
    end if

  end function median


  subroutine sort( vec, len )
    !--------------------------------------------------------------------
    ! This subroutine receives an array vec() and sorts it into ascending
    ! order.
    ! Copied from https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sort.f90
    !--------------------------------------------------------------------
    real, dimension(:), intent(inout)     :: vec
    integer, intent(in)                   :: len
    integer                               :: i
    integer                               :: location

    do i = 1, len-1                           ! eveccept for the last
      location = find_minimum( vec, i, len )  ! find min from this to last
      call swap( vec(i), vec(location) )      ! swap this and the minimum
    end do

  end subroutine sort


  subroutine swap( a, b )
    !--------------------------------------------------------------------
    ! This subroutine swaps the values of its two formal arguments.
    ! Copied from https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sort.f90
    !--------------------------------------------------------------------
    real, intent(inout) :: a, b
    real                :: tmp

    tmp = a
    a   = b
    b   = tmp

  end subroutine  swap


  function find_minimum( vec, start, end ) result( out )
    !--------------------------------------------------------------------
    ! This function returns the location of the minimum in the section
    ! between start and end.
    ! Copied from https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sort.f90
    !--------------------------------------------------------------------
    real, dimension(:), intent(in)     :: vec
    integer, intent(in)                :: start, end
    real                               :: minimum
    integer                            :: location
    integer                            :: i
    integer                            :: out

    minimum  = vec(start)          ! assume the first is the min
    location = start               ! record its position
    do i = start+1, end            ! start with nevect elements
      if (vec(i) < minimum) THEN  ! if vec(i) less than the min?
        minimum  = vec(i)        ! Yes, a new minimum found
        location = i             ! record its position
      end if
    end do
    out = location                 ! return the position

  end function find_minimum


  function findroot_quadratic( a, b, c ) result( root )
    !-----------------------------------------------------------------------
    ! Returns the solution for a quadratic function:
    ! a + bx + cx**2 = 0
    ! Per default returns root2 
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: a, b, c

    ! function return variable
    real, dimension(2) :: root

    ! local variables
    real :: d

    d = b*b - 4.0*a*c
    if (d >= 0.0) then              ! is it solvable?
      d = sqrt(d)
      root(1) = (-b + d)/(2.0*c)    ! first root
      root(2) = (-b - d)/(2.0*c)    ! second root
    else                            ! complex roots
      ! stop 'findroot_quadratic(): There is no real root.'
    end if

  end function findroot_quadratic


  function dgcos( x ) result( dgcos_out )
    !----------------------------------------------------------------   
    ! Calculates the cosine of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: x  ! angle, degrees (0-360)

    ! function return value
    real :: dgcos_out ! cosine value of x when x is in degrees

    !dgcos = dcos(x*pi/180.0)
    dgcos_out = cos(x*pi/180.0)  ! xxx use cos with single-precision compilation

  end function dgcos


  function dgsin( x ) result( dgsin_out )
    !----------------------------------------------------------------   
    ! Calculates the sinus of an angle given in degrees. Equal to 
    ! 'dsin' in Python version.
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: x  ! angle, degrees (0-360)

    ! function return value
    real :: dgsin_out ! sinus value of x when x is in degrees

    !dgsin_out = dsin(x*pi/180.0)
    dgsin_out = sin(x*pi/180.0)   ! xxx use cos with single-precision compilation

  end function dgsin


  function degrees( x ) result( degrees_out )
    !----------------------------------------------------------------   
    ! Returns corresponding degrees if x is given in radians
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: x  ! angle, radians

    ! function return value
    real :: degrees_out

    degrees_out = x*180.0/pi

  end function degrees


  function radians( x ) result( radians_out )
    !----------------------------------------------------------------   
    ! Returns corresponding radians if x is given in degrees
    !----------------------------------------------------------------   
    ! arguments
    real, intent(in) :: x  ! angle, radians

    ! function return value
    real :: radians_out

    radians_out = x*pi/180.0

  end function radians


end module md_sofunutils
