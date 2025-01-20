module md_params_siml_biomee
  !////////////////////////////////////////////////////////////////
  ! Module for handling simulation parameters for BiomeE.
  !----------------------------------------------------------------
  use md_params_core, only: steering_parameters

  implicit none

  private
  public paramstype_siml_biomee

  !----------------------------------------------------------------
  ! Derived type for simulation parameters
  !----------------------------------------------------------------
  type paramstype_siml_biomee

    type(steering_parameters) :: steering
    logical :: do_U_shaped_mortality
    logical :: do_closedN_run
    character(len=30) :: method_photosynth
    character(len=30) :: method_mortality

  end type paramstype_siml_biomee

end module md_params_siml_biomee


