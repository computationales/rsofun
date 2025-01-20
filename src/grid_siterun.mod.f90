module md_grid
  !////////////////////////////////////////////////////////////////
  ! Module for handling variables defining the model grid. For site-
  ! scale simulation this is trivial but still necesssary.
  !----------------------------------------------------------------
  implicit none

  private
  public gridtype

  type gridtype
    integer :: ilon
    integer :: ilat
    real :: lon
    real :: lat
    real :: elv
    real :: nu               ! true anomaly (orbital parameter), recalculated each year for each gridcell in solar()
    real :: lambda           ! true longitude (orbital parameter), recalculated each year for each gridcell in solar()
    real :: decl_angle       ! declination angle (degrees)
    real :: dayl             ! day length (s), is updated daily in waterbal SR
  end type gridtype

end module md_grid
