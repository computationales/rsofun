module pmodel_mod
  !////////////////////////////////////////////////////////////////
  ! Contains main subroutine handling I/O with C and R.
  !----------------------------------------------------------------
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char, c_bool

  implicit none

  private
  public :: pmodel_f

contains

  subroutine pmodel_f(         &
    spinup,                    &   
    spinupyears,               &        
    recycle,                   &    
    firstyeartrend,            &           
    nyeartrend,                &  
    secs_per_tstep,            &     
    in_ppfd,                   &    
    in_netrad,                 &      
    outdt,                     &  
    ltre,                      & 
    ltne,                      & 
    ltrd,                      & 
    ltnd,                      & 
    lgr3,                      & 
    lgn3,                      & 
    lgr4,                      & 
    longitude,                 &      
    latitude,                  &     
    altitude,                  &   
    whc,                       &
    nt,                        &
    par,                       &
    forcing,                   &
    output                     &
    ) bind(C, name = "pmodel_f_")

    !////////////////////////////////////////////////////////////////
    ! Main subroutine to handle I/O with C and R. 
    ! Receives simulation parameters, site parameters, and the full 
    ! simulation's forcing as time series
    !----------------------------------------------------------------
    use md_forcing_pmodel, only: getclimate, getco2, getfapar, get_fpc_grid
    use md_interface_pmodel, only: interfacetype_biosphere, outtype_biosphere, myinterface
    use md_params_core
    use md_biosphere_pmodel, only: biosphere_annual

    implicit none

    ! arguments
    integer(kind=c_int),  intent(in) :: spinup         ! logical type is not supported in the C interface (LTO)
    integer(kind=c_int),  intent(in) :: spinupyears
    integer(kind=c_int),  intent(in) :: recycle
    integer(kind=c_int),  intent(in) :: firstyeartrend
    integer(kind=c_int),  intent(in) :: nyeartrend
    integer(kind=c_int),  intent(in) :: secs_per_tstep
    integer(kind=c_int),  intent(in) :: in_ppfd        ! logical
    integer(kind=c_int),  intent(in) :: in_netrad      ! logical
    integer(kind=c_int),  intent(in) :: outdt
    integer(kind=c_int), intent(in) :: ltre            ! logical
    integer(kind=c_int), intent(in) :: ltne            ! logical
    integer(kind=c_int), intent(in) :: ltrd            ! logical
    integer(kind=c_int), intent(in) :: ltnd            ! logical
    integer(kind=c_int), intent(in) :: lgr3            ! logical
    integer(kind=c_int), intent(in) :: lgn3            ! logical
    integer(kind=c_int), intent(in) :: lgr4            ! logical
    real(kind=c_double),  intent(in) :: longitude
    real(kind=c_double),  intent(in) :: latitude
    real(kind=c_double),  intent(in) :: altitude
    real(kind=c_double),  intent(in) :: whc
    integer(kind=c_int),  intent(in) :: nt ! number of time steps
    real(kind=c_double),  dimension(9), intent(in) :: par  ! free (calibratable) model parameters
    real(kind=c_double),  dimension(nt,12), intent(in) :: forcing  ! array containing all temporally varying forcing data (rows: time steps; columns: 1=air temperature, 2=rainfall, 3=vpd, 4=ppfd, 5=net radiation, 6=sunshine fraction, 7=snowfall, 8=co2, 9=fapar, 10=patm, 11=tmin, 12=tmax) 
    real(kind=c_double),  dimension(nt,19), intent(out) :: output

    ! local variables
    type(outtype_biosphere) :: out_biosphere  ! holds all the output used for calculating the cost or maximum likelihood function 
    integer :: yr, idx_start, idx_end ! npft_local

    !----------------------------------------------------------------
    ! GET SIMULATION PARAMETERS
    !----------------------------------------------------------------
    myinterface%params_siml%steering%do_spinup      = spinup /= 0
    myinterface%params_siml%steering%spinupyears    = spinupyears
    myinterface%params_siml%steering%recycle        = recycle
    myinterface%params_siml%steering%firstyeartrend = firstyeartrend
    myinterface%params_siml%steering%nyeartrend     = nyeartrend

    if (myinterface%params_siml%steering%do_spinup) then
      myinterface%params_siml%steering%runyears = myinterface%params_siml%steering%nyeartrend &
              + myinterface%params_siml%steering%spinupyears
    else
      myinterface%params_siml%steering%runyears = myinterface%params_siml%steering%nyeartrend
      myinterface%params_siml%steering%spinupyears = 0
    endif
    
    myinterface%params_siml%in_ppfd            = in_ppfd /= 0
    myinterface%params_siml%in_netrad          = in_netrad /= 0
    myinterface%params_siml%outdt              = outdt
    myinterface%params_siml%ltre               = ltre /= 0
    myinterface%params_siml%ltne               = ltne /= 0
    myinterface%params_siml%ltrd               = ltrd /= 0
    myinterface%params_siml%ltnd               = ltnd /= 0
    myinterface%params_siml%lgr3               = lgr3 /= 0
    myinterface%params_siml%lgn3               = lgn3 /= 0
    myinterface%params_siml%lgr4               = lgr4 /= 0
    myinterface%params_siml%secs_per_tstep     = secs_per_tstep

    !----------------------------------------------------------------
    ! GET GRID INFORMATION
    !----------------------------------------------------------------
    myinterface%grid%lon = real( longitude )
    myinterface%grid%lat = real( latitude )
    myinterface%grid%elv = real( altitude )

    !----------------------------------------------------------------
    ! GET SOIL PARAMETERS
    !----------------------------------------------------------------
    myinterface%whc_prescr = real( whc )
    
    !----------------------------------------------------------------
    ! GET CALIBRATABLE MODEL PARAMETERS (so far a small list)
    !----------------------------------------------------------------
    myinterface%params_calib%kphio              = real(par(1))
    myinterface%params_calib%kphio_par_a        = real(par(2))
    myinterface%params_calib%kphio_par_b        = real(par(3))
    myinterface%params_calib%soilm_thetastar    = real(par(4))
    myinterface%params_calib%soilm_betao        = real(par(5))
    myinterface%params_calib%beta_unitcostratio = real(par(6))
    myinterface%params_calib%rd_to_vcmax        = real(par(7))
    myinterface%params_calib%tau_acclim         = real(par(8))
    myinterface%params_calib%kc_jmax            = real(par(9))

    !----------------------------------------------------------------
    ! GET VEGETATION COVER (fractional projective cover by PFT)
    !----------------------------------------------------------------
    myinterface%fpc_grid(:) = get_fpc_grid( myinterface%params_siml )
    
    do yr=1,myinterface%params_siml%steering%runyears

      !----------------------------------------------------------------
      ! Define simulations "steering" variables (forcingyear, etc.)
      !----------------------------------------------------------------
      myinterface%steering = get_steering( yr, myinterface%params_siml%steering )

      !----------------------------------------------------------------
      ! Get external (environmental) forcing
      !----------------------------------------------------------------
      ! Get climate variables for this year (full fields and 365 daily values for each variable)
      myinterface%climate(:) = getclimate(nt, &
                                          forcing, &
                                          myinterface%steering%climateyear_idx, &
                                          myinterface%params_siml%in_ppfd,  &
                                          myinterface%params_siml%in_netrad &
                                          )

      ! Get annual, gobally uniform CO2
      myinterface%pco2 = getco2(  nt, &
                                  forcing, &
                                  myinterface%steering%forcingyear, &
                                  myinterface%params_siml%steering%firstyeartrend &
                                  )

      !----------------------------------------------------------------
      ! Get prescribed fAPAR if required (otherwise set to dummy value)
      !----------------------------------------------------------------
      myinterface%vegcover(:) = getfapar( &
                                        nt, &
                                        forcing, &
                                        myinterface%steering%forcingyear_idx &
                                        )

      !----------------------------------------------------------------
      ! Call biosphere (wrapper for all modules, contains gridcell loop)
      !----------------------------------------------------------------
      out_biosphere = biosphere_annual() 
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      ! Populate Fortran output array which is passed back to C/R
      !----------------------------------------------------------------
      if (yr > myinterface%params_siml%steering%spinupyears ) then

        idx_start = (myinterface%steering%forcingyear_idx - 1) * ndayyear + 1
        idx_end   = idx_start + ndayyear - 1

        output(idx_start:idx_end,1)  = dble(out_biosphere%fapar(:))  
        output(idx_start:idx_end,2)  = dble(out_biosphere%gpp(:))    
        output(idx_start:idx_end,3)  = dble(out_biosphere%transp(:)) 
        output(idx_start:idx_end,4)  = dble(out_biosphere%latenth(:))
        output(idx_start:idx_end,5)  = dble(out_biosphere%pet(:))
        output(idx_start:idx_end,6)  = dble(out_biosphere%vcmax(:))  
        output(idx_start:idx_end,7)  = dble(out_biosphere%jmax(:))    
        output(idx_start:idx_end,8)  = dble(out_biosphere%vcmax25(:)) 
        output(idx_start:idx_end,9)  = dble(out_biosphere%jmax25(:))
        output(idx_start:idx_end,10) = dble(out_biosphere%gs_accl(:))
        output(idx_start:idx_end,11) = dble(out_biosphere%wscal(:))
        output(idx_start:idx_end,12) = dble(out_biosphere%chi(:))
        output(idx_start:idx_end,13) = dble(out_biosphere%iwue(:))
        output(idx_start:idx_end,14) = dble(out_biosphere%rd(:))
        output(idx_start:idx_end,15) = dble(out_biosphere%tsoil(:))
        output(idx_start:idx_end,16) = dble(out_biosphere%netrad(:))
        output(idx_start:idx_end,17) = dble(out_biosphere%wcont(:))
        output(idx_start:idx_end,18) = dble(out_biosphere%snow(:))
        output(idx_start:idx_end,19) = dble(out_biosphere%cond(:))

      end if

    enddo

  end subroutine pmodel_f

end module pmodel_mod
