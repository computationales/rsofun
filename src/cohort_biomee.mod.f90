module md_cohort
  !////////////////////////////////////////////////////////////////
  ! Module containing BiomeE cohort definitions
  !----------------------------------------------------------------
  use md_interface_biomee, only: myinterface, spec_data_type, MAX_LEVELS
  use md_params_core, only: pi
  use md_classdefs

  ! define data types and constants
  implicit none
  private
  !=============== Public types ===========================================================
  public :: cohort_type

  !=============== Public subroutines =====================================================


  integer, public, parameter :: LEAF_OFF                 = 0
  integer, public, parameter :: LEAF_ON                  = 1
  integer, public, parameter :: PT_C3                    = 0                     ! physiology types
  integer, public, parameter :: PT_C4                    = 1                     ! physiology types

  !=============== Cohort level data type =============================================================
  type :: cohort_type

    !======================== State
    ! The cohort state contains cohort properties which persist across the years

    !===== Metadata
    integer :: ccID       = -1           ! cohort ID
    integer :: layer      = 1            ! the layer of this cohort (numbered from top, top layer=1)
    integer :: firstlayer = 0            ! 0 = never been in the first layer; 1 = at least one year in first layer
    integer :: species    = 1            ! vegetation species

    !===== Population structure
    real :: nindivs   = 1.0              ! density of vegetation, tree/m2
    real :: age       = 0.0              ! age of cohort, years
    real :: topyear   = 0.0              ! number of years the cohort is in top layer
    real :: dbh       = 0.0              ! diameter at breast height, m
    real :: height    = 0.0              ! vegetation height, m
    real :: crownarea = 1.0              ! crown area, m2 tree-1
    real :: leafarea  = 0.0              ! total area of leaves, m2 tree-1

    !===== Biological prognostic variables
    real    :: gdd        = 0.0          ! growing degree-day (phenology)
    integer :: status     = LEAF_OFF     ! growth status of plant
    real    :: leaf_age   = 0.0          ! leaf age (years)

    !===== Organic pools
    type(orgpool) :: pleaf               ! leaf biomass, kg tree-1
    type(orgpool) :: proot               ! root biomass, kg tree-1
    type(orgpool) :: psapw               ! sapwood biomass, kg tree-1
    type(orgpool) :: pwood               ! heartwood (non-living) biomass, kg tree-1
    type(orgpool) :: pseed               ! biomass put aside for future progeny, kg tree-1
    type(orgpool) :: plabl               ! labile pool, temporary storage of N and C, kg tree-1

    !===== Nitrogen model related variables (persistent from year to year)
    real    :: bl_max             = 0.0           ! Max. leaf biomass, kg C tree-1
    real    :: br_max             = 0.0           ! Max. fine root biomass, kg C tree-1

    !=================== Temporary variables
    ! Contrary to the state vairbales, temporary variables are reset every step, day, or year (as appropriate).

    !===== Fast step fluxes, kg timestep-1 tree-1
    type(common_fluxes) :: fast_fluxes
    real    :: resl              = 0.0            ! leaf respiration, kg C timestep-1 tree-1
    real    :: resr              = 0.0            ! root respiration, kg C timestep-1 tree-1

    !===== Daily fluxes, kg day-1 tree-1
    type(common_fluxes) :: daily_fluxes

    !===== Annual fluxes, kg yr-1 tree-1
    type(common_fluxes) :: annual_fluxes
    real    :: NPPleaf            = 0.0           ! C allocated to leaf, kg C yr-1 tree-1
    real    :: NPProot            = 0.0           ! C allocated to root, kg C yr-1 tree-1
    real    :: NPPwood            = 0.0           ! C allocated to wood, kg C yr-1 tree-1

    !===== Annual fluxes due to tree death, kg m-2 yr-1
    real    :: n_deadtrees        = 0.0           ! plant to soil N flux due to mortality (kg N m-2 yr-1)
    real    :: c_deadtrees        = 0.0           ! plant to soil C flux due to mortality (kg C m-2 yr-1)
    real    :: m_turnover         = 0.0           ! C turnover due to mortality and tissue turnover (kg C m-2 yr-1)
    real    :: deathrate          = 0.0           ! Deathrate (0 to 1)

    !===== Water uptake-related variables
    real    :: WupL(MAX_LEVELS)   = 0.0           ! normalized vertical distribution of uptake

    !===== Photosynthesis variables
    real    :: An_op              = 0.0           ! mol C/(m2 of leaf per year)
    real    :: An_cl              = 0.0           ! mol C/(m2 of leaf per year)
    real    :: C_growth           = 0.0           ! Carbon gain since last growth, kg C day-1 tree-1
    real    :: N_growth           = 0.0           ! Nitrogen used for plant tissue growth, kg N day-1 tree-1
    real    :: resg               = 0.0           ! growth respiration, kg C day-1 tree-1

    !===== Memory variables used for computing deltas
    real    :: DBH_ys            = 0.0            ! DBH at the begining of a year (growing season)
    real    :: BA_ys             = 0.0            ! Basal area at the beginning og a year
    
    contains

      !========= Derived variables
      ! Variables which are function of any state or temporary variable defined above.
      ! They can be used like any normal variable except that they are followed by parenthesis
      ! (indicating that they are being evaluated every time).
      ! The advantage is that they cannot be out-of-sync with the underlying data, at the cost of a negligeable
      ! increase in the number of operations.

      procedure NSNmax
      procedure W_supply
      procedure rootareaL
      procedure lai
      procedure basal_area
      procedure volume
      procedure merge_in
      procedure layerfrac
      procedure sp

      !========== Other member procedures

      procedure reset_cohort

  end type cohort_type

contains

  function sp(self) result(res)
    type(spec_data_type) :: res
    class(cohort_type) :: self

    res = myinterface%params_species(self%species)
  end function sp

  subroutine merge_in(self, other)
    !////////////////////////////////////////////////////////////////
    ! Merge other cohort into self
    !---------------------------------------------------------------
    class(cohort_type), intent(inout) :: self
    type(cohort_type), intent(in) :: other
    real :: x1, x2 ! normalized relative weights

    if (other%nindivs > 0.0 .or. self%nindivs > 0.0) then

      x1 = other%nindivs / (other%nindivs + self%nindivs)
      x2 = self%nindivs / (other%nindivs + self%nindivs)

      ! update number of individuals in merged cohort
      self%nindivs = other%nindivs + self%nindivs

      ! Carbon
      self%pleaf%c%c12 = x1*other%pleaf%c%c12 + x2*self%pleaf%c%c12
      self%proot%c%c12 = x1*other%proot%c%c12 + x2*self%proot%c%c12
      self%psapw%c%c12 = x1*other%psapw%c%c12 + x2*self%psapw%c%c12
      self%pwood%c%c12 = x1*other%pwood%c%c12 + x2*self%pwood%c%c12
      self%pseed%c%c12 = x1*other%pseed%c%c12 + x2*self%pseed%c%c12
      self%plabl%c%c12 = x1*other%plabl%c%c12 + x2*self%plabl%c%c12

      ! Nitrogen
      self%pleaf%n%n14 = x1*other%pleaf%n%n14 + x2*self%pleaf%n%n14
      self%proot%n%n14 = x1*other%proot%n%n14 + x2*self%proot%n%n14
      self%psapw%n%n14 = x1*other%psapw%n%n14 + x2*self%psapw%n%n14
      self%pwood%n%n14 = x1*other%pwood%n%n14 + x2*self%pwood%n%n14
      self%pseed%n%n14 = x1*other%pseed%n%n14 + x2*self%pseed%n%n14
      self%plabl%n%n14 = x1*other%plabl%n%n14 + x2*self%plabl%n%n14

      ! Cohort age
      self%age = x1*other%age + x2*self%age
      self%topyear = x1*other%topyear + x2*self%topyear
    endif

  end subroutine merge_in

  subroutine reset_cohort(self)
    ! Reset cohort temporary data (used yearly)
    class(cohort_type), intent(inout) :: self

    ! Save last year's values
    self%DBH_ys        = self%dbh
    self%BA_ys         = basal_area(self)

    self%WupL(:)       = 0.0

    self%An_op         = 0.0
    self%An_cl         = 0.0
    self%N_growth      = 0.0
    self%N_growth      = 0.0

    self%resl          = 0.0
    self%resr          = 0.0
    self%resg          = 0.0

    self%fast_fluxes   = common_fluxes()

    ! Reset daily
    self%daily_fluxes  = common_fluxes()

    ! annual
    self%annual_fluxes = common_fluxes()
    self%NPPleaf       = 0.0
    self%NPProot       = 0.0
    self%NPPwood       = 0.0

    self%n_deadtrees   = 0.0
    self%c_deadtrees   = 0.0
    self%m_turnover    = 0.0
    self%deathrate     = 0.0
  end subroutine reset_cohort

  function NSNmax(self) result(res)
    real :: res
    class(cohort_type) :: self

    ! Local variable
    type(spec_data_type) :: sp

    sp = self%sp()

    res = sp%fNSNmax * &
            (self%bl_max / (sp%CNleaf0 * sp%leafLS) + self%br_max / sp%CNroot0)
  end function NSNmax

  function W_supply(self) result(res)
    ! potential water uptake rate per unit time per tree
    real :: res
    class(cohort_type) :: self

    res = sum(self%WupL(:))
  end function W_supply

  function rootareaL(self, level) result(res)
    ! Root length per layer, m of root/m
    real :: res
    integer :: level
    class(cohort_type) :: self

    ! Local variable
    type(spec_data_type) :: sp

    sp = self%sp()

    res = rootarea(self) * sp%root_frac(level)
  end function rootareaL

  function rootarea(self) result(res)
    ! total fine root area per tree
    real :: res
    class(cohort_type) :: self

    ! Local variable
    type(spec_data_type) :: sp

    sp = self%sp()

    res  = self%proot%c%c12 * sp%SRA
  end function rootarea

  function lai(self) result(res)
    ! Leaf area index: surface of leaves per m2 of crown
    real :: res
    class(cohort_type) :: self

    res = self%leafarea / self%crownarea !(cohort%crownarea *(1.0-sp%internal_gap_frac))
  end function lai

  function basal_area(self) result(res)
    ! Tree basal area, m2 tree-1
    real :: res
    class(cohort_type) :: self

    res = pi/4 * self%dbh * self%dbh
  end function basal_area

  function volume(self) result(res)
    ! Tree basal volume, m3 tree-1
    real :: res
    class(cohort_type) :: self
    ! Local variable
    type(spec_data_type) :: sp

    sp = self%sp()

    res = (self%psapw%c%c12 + self%pwood%c%c12) / sp%rho_wood
  end function volume

  function layerfrac(self) result(res)
    ! Fraction of layer area occupied by this cohort
    real :: res
    class(cohort_type) :: self

    res = self%nindivs * self%crownarea
  end function layerfrac

end module md_cohort
