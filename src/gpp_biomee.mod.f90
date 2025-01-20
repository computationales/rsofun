module md_gpp_biomee
  !////////////////////////////////////////////////////////////////
  ! Module containing a wrapper for using different photosynthesis
  ! schemes within BiomeE.
  ! Code for gs_leuning photosynthesis option is adopted from BiomeE https://doi.org/10.5281/zenodo.7125963.
  ! Code for pmodel photosynthesis option is for P-model (Stocker et al., 2020 GMD)
  !----------------------------------------------------------------
  use datatypes_biomee
  use md_interface_biomee, only: myinterface
  use md_soil_biomee, only: water_supply_layer
  use md_sofunutils, only: calc_esat

  implicit none

  private
  public gpp, getpar_modl_gpp

  !-----------------------------------------------------------------------
  ! P-model parameters created here for pmodel option. takes no effect in gs_leuning option
  !-----------------------------------------------------------------------
  type paramstype_gpp
    real :: beta         ! Unit cost of carboxylation (dimensionless)
    real :: soilm_thetastar
    real :: soilm_betao
    real :: rd_to_vcmax  ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
    real :: tau_acclim   ! acclimation time scale of photosynthesis (d)
    real :: kc_jmax

    ! these should be species-specific, temporary solution to put them here
    real :: kphio        ! quantum yield efficiency at optimal temperature, phi_0 (Stocker et al., 2020 GMD Eq. 10)
    real :: kphio_par_a  ! shape parameter of temperature-dependency of quantum yield efficiency
    real :: kphio_par_b  ! optimal temperature of quantum yield efficiency (deg C)

  end type paramstype_gpp

  type(paramstype_gpp) :: params_gpp

contains

  subroutine gpp( forcing, vegn )
    !//////////////////////////////////////////////////////////////////////
    ! GPP
    ! Calculates light availability and photosynthesis for each cohort 
    ! and sets the following cohort-level variables:
    ! - An_op   
    ! - An_cl   
    ! - w_scale 
    ! - transp  
    ! - gpp     
    !
    ! Subroutines from BiomeE-Allocation
    !------------------------------------------------------------------------
    use md_forcing_biomee, only: climate_type
    use md_photosynth, only: pmodel, zero_pmodel, outtype_pmodel, calc_ftemp_inst_rd
    use md_photosynth, only: calc_kphio_temp, calc_soilmstress
    use md_params_core, only: kTkelvin, kfFEC, c_molmass
    use md_sofunutils, only: dampen_variability

    type(climate_type), intent(in):: forcing
    type(vegn_tile_type), intent(inout) :: vegn

    ! local variables used for BiomeE-Allocation part
    type(cohort_type), pointer :: cc
    type(cohort_item), pointer :: it
    integer :: i
    real   :: rad_top                                      ! downward radiation at the top of the canopy, W/m2
    real   :: rad_net                                      ! net radiation absorbed by the canopy, W/m2
    real   :: Tair, TairK                                  ! air temperature, degC and degK
    real   :: cana_q                                       ! specific humidity in canopy air space, kg/kg
    real   :: cana_co2                                     ! co2 concentration in canopy air space, mol CO2/mol dry air
    real   :: p_surf                                       ! surface pressure, Pa
    real   :: water_supply                                 ! water supply per m2 of leaves
    real   :: fw, fs                                       ! wet and snow-covered fraction of leaves
    real   :: psyn                                         ! net photosynthesis, mol C/(m2 of leaves s)
    real   :: resp                                         ! leaf respiration, mol C/(m2 of leaves s)
    real   :: w_scale2, transp                             ! mol H20 per m2 of leaf per second
    real   :: f_light(nlayers_max+1)                       ! incident light fraction at top of a given layer
    real   :: LAIlayer(nlayers_max)                        ! leaf area index per layer, corrected for gaps (representative for the tree-covered fraction)
    real   :: accuCAI
    real   :: par                                          ! just for temporary use
    real, allocatable :: fapar_tree(:)                     ! tree-level fAPAR based on LAI within the crown
    real, dimension(nlayers_max-1) :: fapar_layer
    real, parameter :: kappa = 0.5                         ! light extinction coefficient of crown layers
    real, parameter :: f_gap = 0.1

    ! local variables used for P-model part
    real :: kphio_temp
    type(outtype_pmodel) :: out_pmodel      ! list of P-model output variables

    allocate(fapar_tree(vegn%n_cohorts())) ! no need to deallocate
    fapar_layer(:) = 0.0
    LAIlayer(:) = 0.0

    !-----------------------------------------------------------
    ! Canopy light absorption
    !-----------------------------------------------------------
    ! ! Calculate kappa according to sun zenith angle 
    ! kappa = extinct/max(cosz,0.01)

    ! Sum leaf area over cohorts in each crown layer -> LAIlayer(layer)
    i = 0
    it => vegn%heap
    do while (associated(it))
      cc => it%cohort
      i = i + 1
      LAIlayer(cc%layer) = LAIlayer(cc%layer) + cc%leafarea * cc%nindivs / (1.0 - f_gap)
      fapar_tree(i) = 1.0 - exp(-kappa * cc%leafarea / cc%crownarea)   ! at individual-level: cc%leafarea represents leaf area index within the crown
      fapar_layer(cc%layer) = fapar_layer(cc%layer) + fapar_tree(i) * cc%crownarea * cc%nindivs
      it => it%next
    end do

    ! Get light received at each crown layer as a fraction of top-of-canopy -> f_light(layer) 
    f_light(:) = 0.0
    f_light(1) = 1.0
    do i = 2, nlayers_max
      ! f_light(i) = f_light(i-1) * (exp(-kappa * LAIlayer(i-1)) + f_gap)                     ! originally in LM3-PPA
      ! f_light(i) = f_light(i-1) * ((1.0 - f_gap) * exp(-kappa * LAIlayer(i-1)) + f_gap)     ! corrected version, corresponding to original LM3-PPA approach
      f_light(i) = f_light(i-1) * (1.0 - fapar_layer(i-1))                                    ! alternative version for conserving energy
    enddo


    if (trim(myinterface%params_siml%method_photosynth) == "gs_leuning") then
      !===========================================================
      ! Original BiomeE-Allocation
      !-----------------------------------------------------------
      ! Water supply for photosynthesis, Layers
      call water_supply_layer(vegn)

      ! Photosynthesis
      accuCAI = 0.0

      it => vegn%heap
      do while (associated(it))

        cc => it%cohort
        associate ( sp => cc%sp() )

        if (cc%status == LEAF_ON) then   !.and. cc%lai > 0.1

          ! Convert forcing data
          rad_top  = f_light(cc%layer) * forcing%radiation ! downward radiation at the top of the canopy, W/m2

          !===============================
          ! ORIGINAL
          !===============================
          rad_net  = f_light(cc%layer) * forcing%radiation * 0.9 ! net radiation absorbed by the canopy, W/m2
          p_surf   = forcing%P_air  ! Pa
          TairK    = forcing%Tair ! K
          Tair     = forcing%Tair - 273.16 ! degC
          cana_q   = (calc_esat(Tair) * forcing%RH * h2o_molmass) / (p_surf * kMa)  ! air specific humidity, kg/kg
          cana_co2 = forcing%CO2 ! co2 concentration in canopy air space, mol CO2/mol dry air

          ! recalculate the water supply to mol H20 per m2 of leaf per second
          water_supply = cc%W_supply() / (cc%leafarea * myinterface%step_seconds * h2o_molmass * 1e-3) ! mol m-2 leafarea s-1

          !call get_vegn_wet_frac (cohort, fw=fw, fs=fs)
          fw = 0.0
          fs = 0.0

          call gs_leuning(rad_top, rad_net, TairK, cana_q, cc%lai(), &
            p_surf, water_supply, cc%species, sp%pt, &
            cana_co2, extinct, fs+fw, &
            psyn, resp, w_scale2, transp )

          ! store the calculated photosynthesis, photorespiration, and transpiration for future use in growth
          cc%An_op   = psyn   ! molC s-1 m-2 of leaves ! net photosynthesis, mol C/(m2 of leaves s)
          cc%An_cl   = -resp  ! molC s-1 m-2 of leaves
          cc%fast_fluxes%trsp = transp * h2o_molmass * 1e-3 * cc%leafarea * myinterface%step_seconds      ! Transpiration (kgH2O/(tree step), Weng, 2017-10-16
          cc%resl = -resp * c_molmass * 1e-3 * cc%leafarea * myinterface%step_seconds ! kgC tree-1 step-1
          cc%fast_fluxes%gpp = (psyn - resp) * c_molmass * 1e-3 * cc%leafarea * myinterface%step_seconds ! kgC step-1 tree-1

          else

          ! no leaves means no photosynthesis and no stomatal conductance either
            cc%An_op   = 0.0
            cc%An_cl   = 0.0
            cc%fast_fluxes%gpp   = 0.0
            cc%fast_fluxes%trsp  = 0.0

          endif
        end associate

        it => it%next
      end do

    else if (trim(myinterface%params_siml%method_photosynth) == "pmodel") then
      !===========================================================
      ! P-model
      !-----------------------------------------------------------
      ! Calculate environmental conditions with memory, time scale 
      ! relevant for Rubisco turnover
      !----------------------------------------------------------------
      if (.not. vegn%dampended_forcing%initialized) then
        vegn%dampended_forcing%co2  = forcing%CO2 * 1.0e6
        vegn%dampended_forcing%temp = (forcing%Tair - kTkelvin)
        vegn%dampended_forcing%vpd  = forcing%vpd
        vegn%dampended_forcing%patm = forcing%P_air
        vegn%dampended_forcing%par = forcing%radiation
        vegn%dampended_forcing%initialized = .True.
      else
        vegn%dampended_forcing%co2  = dampen_variability( forcing%CO2 * 1.0e6,        params_gpp%tau_acclim, &
                vegn%dampended_forcing%co2 )
        vegn%dampended_forcing%temp = dampen_variability( (forcing%Tair - kTkelvin),  params_gpp%tau_acclim, &
                vegn%dampended_forcing%temp)
        vegn%dampended_forcing%vpd  = dampen_variability( forcing%vpd,                &
                params_gpp%tau_acclim, vegn%dampended_forcing%vpd )
        vegn%dampended_forcing%patm = dampen_variability( forcing%P_air,              &
                params_gpp%tau_acclim, vegn%dampended_forcing%patm )
        vegn%dampended_forcing%par = dampen_variability( forcing%radiation,              &
                params_gpp%tau_acclim, vegn%dampended_forcing%par )
      end if

      !----------------------------------------------------------------
      ! Instantaneous temperature effect on quantum yield efficiency
      !----------------------------------------------------------------
      kphio_temp = calc_kphio_temp( (forcing%Tair - kTkelvin), &
              .false., &    ! no C4
              params_gpp%kphio, &
              params_gpp%kphio_par_a, &
              params_gpp%kphio_par_b )

      !----------------------------------------------------------------
      ! Photosynthesis for each cohort
      !----------------------------------------------------------------
      i = 0
      it => vegn%heap
      do while (associated(it))

        cc => it%cohort
        i = i + 1

        if (cc%status == LEAF_ON) then

          ! photosynthetically active radiation level at this layer
          par = f_light(cc%layer) * vegn%dampended_forcing%par * kfFEC * 1.0e-6

          !----------------------------------------------------------------
          ! P-model call for C3 plants to get a list of variables that are 
          ! acclimated to slowly varying conditions
          !----------------------------------------------------------------
          out_pmodel = pmodel(  &
                                kphio          = kphio_temp, &
                                beta           = params_gpp%beta, &
                                kc_jmax        = params_gpp%kc_jmax, &
                                ppfd           = par, &
                                co2            = vegn%dampended_forcing%co2, &
                                tc             = vegn%dampended_forcing%temp, &
                                vpd            = vegn%dampended_forcing%vpd, &
                                patm           = vegn%dampended_forcing%patm, &
                                c4             = .false., &
                                method_optci   = "prentice14", &
                                method_jmaxlim = "wang17" &
                                )

          ! irrelevant variables for this setup  
          cc%An_op   = 0.0
          cc%An_cl   = 0.0
          cc%fast_fluxes%trsp  = 0.0

          ! quantities per tree and cumulated over seconds in time step (kgC step-1 tree-1 )
          cc%fast_fluxes%gpp = par * fapar_tree(i) * out_pmodel%lue * cc%crownarea * myinterface%step_seconds * 1.0e-3
          cc%resl = fapar_tree(i) * out_pmodel%vcmax25 * params_gpp%rd_to_vcmax * calc_ftemp_inst_rd( forcing%Tair - kTkelvin ) &
            * cc%crownarea * myinterface%step_seconds * c_molmass * 1.0e-3

        else

          ! no leaves means no photosynthesis and no stomatal conductance either
          cc%An_op   = 0.0
          cc%An_cl   = 0.0
          cc%resl    = 0.0
          cc%fast_fluxes%gpp   = 0.0
          cc%fast_fluxes%trsp  = 0.0

        endif

        it => it%next
      end do

    end if

  end subroutine gpp


  subroutine gs_leuning( rad_top, rad_net, tl, ea, lai, &
    p_surf, ws, species, pt, ca, kappa, leaf_wet, &
    apot, acl,w_scale2, transp )

    ! taking from params_core (SOFUN)
    use md_params_core, only: kR, kTkelvin
    use, intrinsic :: ieee_arithmetic

    real,    intent(in)    :: rad_top ! PAR dn on top of the canopy, w/m2
    real,    intent(in)    :: rad_net ! PAR net on top of the canopy, w/m2
    real,    intent(in)    :: tl   ! leaf temperature, degK
    real,    intent(in)    :: ea   ! specific humidity in the canopy air, kg/kg
    real,    intent(in)    :: lai  ! leaf area index
    !real,    intent(in)    :: leaf_age ! age of leaf since budburst (deciduos), days
    real,    intent(in)    :: p_surf ! surface pressure, Pa
    real,    intent(in)    :: ws   ! water supply, mol H20/(m2 of leaf s)
    integer, intent(in)    :: species  ! species
    integer, intent(in)    :: pt   ! physiology type (C3 or C4)
    real,    intent(in)    :: ca   ! concentartion of CO2 in the canopy air space, mol CO2/mol dry air
    real,    intent(in)    :: kappa ! canopy extinction coefficient (move inside f(species))
    real,    intent(in)    :: leaf_wet ! fraction of leaf that's wet or snow-covered
    ! integer, intent(in)    :: layer  ! the layer of this canopy
    ! note that the output is per area of leaf; to get the quantities per area of
    ! land, multiply them by LAI
    !real,    intent(out)   :: gs   ! stomatal conductance, m/s
    real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
    real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
    real,    intent(out)   :: w_scale2, transp  ! transpiration, mol H20/(m2 of leaf s)
    
    ! local variables     
    ! photosynthesis
    real :: vm
    real :: kc, ko ! Michaelis-Menten constants for CO2 and O2, respectively
    real :: ci
    real :: capgam ! CO2 compensation point
    real :: f2, f3
    real :: coef0, coef1
    real :: Resp

    ! conductance related
    real :: gs
    real :: b
    real :: ds  ! humidity deficit, kg/kg
    real :: hl  ! saturated specific humidity at the leaf temperature, kg/kg
    real :: do1

    ! misceleneous
    real :: dum2
    real, parameter :: light_crit = 0
    real, parameter :: gs_lim = 0.25

    ! new average computations
    real :: lai_eq
    real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta 
    real :: light_top
    real :: par_net
    real :: Ag
    real :: An
    real :: Ag_l
    real :: Ag_rb
    real :: anbar
    real :: gsbar
    real :: w_scale
    real :: kk

    ! soil water stress
    real :: Ed, an_w, gs_w


    associate (spdata => myinterface%params_species)

    b = 0.01
    do1 = 0.15 ! kg/kg

    ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
    ! empirical relationship from McCree is light=rn*0.0000046
    light_top = rad_top*rad_phot;
    par_net   = rad_net*rad_phot;

    ! calculate humidity deficit, kg/kg
    hl = qscomp(tl, p_surf)
    ds = max(hl - ea,0.0)

    !  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*kPo/p_surf;
    !  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*kPo/p_surf;
    !  vm=spdata(species)%Vmax*exp(3000.0*(1.0/288.2-1.0/tl));
    ! corrected by Weng, 2013-01-17
    ! Weng, 2013-01-10
    ko=0.248    * exp(35948/kR*(1.0/298.2-1.0/tl))*kPo/p_surf ! Weng, 2013-01-10
    kc=0.000404 * exp(59356/kR*(1.0/298.2-1.0/tl))*kPo/p_surf ! Weng, 2013-01-10
    vm=spdata(species)%Vmax*exp(24920/kR*(1.0/298.2-1.0/tl)) ! / ((layer-1)*1.0+1.0) ! Ea = 33920
    !decrease Vmax due to aging of temperate deciduous leaves 
    !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
    !! Turned off by Weng, 2013-02-01, since we can't trace new leaves
    !  if (spdata(species)%leaf_age_tau>0 .and. leaf_age>spdata(species)%leaf_age_onset) then
    !     vm=vm*exp(-(leaf_age-spdata(species)%leaf_age_onset)/spdata(species)%leaf_age_tau)
    !  endif

    ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986
    capgam=0.5*kc/ko*0.21*0.209; ! Farquhar & Caemmerer 1982

    ! Find respiration for the whole canopy layer

    !  Resp=spdata(species)%gamma_resp*vm*lai /((layer-1)*1.0+1.0)  ! Weng, 2013-01-17 add '/ ((layer-1)*1.0+1.0)'

    ! 2014-09-03, for Nitrogen model: resp = D*(A + B*LMA)
    ! (A+B*LMA) = LNA, D=Vmax/LNA = 25E-6/0.0012 = 0.02 for a standard deciduous species
    !! Leaf resp as a function of nitrogen
    !  Resp=spdata(species)%gamma_resp*0.04*spdata(species)%LNA  & ! basal rate, mol m-2 s-1
    !       * exp(24920/kR*(1.0/298.2-1.0/tl))         & ! temperature scaled
    !       * lai                                        & ! whole canopy
    !       /((layer-1)*1.0+1.0)                         !
    !! as a function of LMA
    !  Resp=(spdata(species)%gamma_LNbase*spdata(species)%LNbase+spdata(species)%gamma_LMA*spdata(species)%LMA)  & ! basal rate, mol m-2 s-1
    !  Resp=spdata(species)%gamma_LNbase*(2.5*spdata(species)%LNA-1.5*spdata(species)%LNbase)     & ! basal rate, mol m-2 s-1
    Resp = spdata(species)%gamma_LN/seconds_per_year & ! per seconds, ! basal rate, mol m-2 s-1
            * spdata(species)%LNA * lai / (c_molmass * 1e-3)    &     ! whole canopy, ! basal rate, mol m-2 s-1
            * exp(24920/kR*(1.0/298.2-1.0/tl))     ! temperature scaled
    !                                  &
    !       /((layer-1)*1.0+1.0)
    ! Temperature effects
    Resp=Resp/((1.0+exp(0.4*(5.0 - tl + kTkelvin)))*(1.0+exp(0.4*(tl - 45.0 - kTkelvin))));


    ! ignore the difference in concentrations of CO2 near
    !  the leaf and in the canopy air, rb=0.
    Ag_l=0.;
    Ag_rb=0.;
    Ag=0.;
    anbar=-Resp/lai;
    gsbar=b;
    ! find the LAI level at which gross photosynthesis rates are equal
    ! only if PAR is positive
    if ( light_top > light_crit ) then

      coef0=(1+ds/do1)/spdata(species)%m_cond
      ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0)

      if (ci>capgam) then

        if (pt==PT_C4) then ! C4 species

          f2=vm
          f3=18000.0*vm*ci ! 18000 or 1800?
          kk = 1.0

        else ! C3 species

          coef1=kc*(1.0+0.209/ko)
          f2=vm*(ci-capgam)/(ci+coef1)
          f3=vm/2.
          kk = (ci+2.*capgam)/(ci-capgam)

        end if

        dum2=min(f2,f3)

        ! find LAI level at which rubisco limited rate is equal to light limited rate
        lai_eq = -log(dum2*kk/(kappa*spdata(species)%alpha_phot*light_top))/kappa
        lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

        ! gross photosynthesis for light-limited part of the canopy
        Ag_l   = spdata(species)%alpha_phot * par_net &
                * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))

        ! We test if Ag_l is NaN and, if so, we set it to 0
        ! (which is the value Ag_l tends to just before exploding due to numerical instability)
        if (ieee_is_nan(Ag_l)) Ag_l = 0.0

        ! gross photosynthesis for rubisco-limited part of the canopy
        Ag_rb  = dum2 * lai_eq

        Ag= (Ag_l + Ag_rb)/((1.0+exp(0.4*(5.0-tl+kTkelvin)))*(1.0+exp(0.4*(tl-45.0-kTkelvin))))
        An=Ag-Resp
        anbar=An/lai

        if (anbar>0.0) then
          gsbar=anbar/(ci-capgam)/coef0
        endif

      endif ! ci>capgam

    endif ! light is available for photosynthesis

    an_w=anbar

    if (an_w > 0.) then
      an_w=an_w*(1-spdata(species)%wet_leaf_dreg*leaf_wet);
    endif
    gs_w = 1.56 * gsbar *(1-spdata(species)%wet_leaf_dreg*leaf_wet); !Weng: 1.56 for H2O?

    if (gs_w > gs_lim) then
      if (an_w > 0.) an_w = an_w*gs_lim/gs_w
      gs_w = gs_lim
    endif

    ! find water availability diagnostic demand
    Ed = gs_w * ds * kMa / h2o_molmass ! ds * kMa * h2o_molmass is the humidity deficit in [mol_h2o/mol_air]

    ! the factor kMa/mol_h2o makes units of gs_w and humidity deficit ds compatible:
    if (Ed>ws) then
      w_scale = ws/Ed
      gs_w = w_scale * gs_w
      if (an_w > 0.0) an_w = an_w * w_scale
      if (an_w < 0.0 .and. gs_w >b) gs_w = b
    endif

    gs=gs_w
    apot=an_w
    acl=-Resp/lai
    transp = min(ws,Ed) ! mol H20/(m2 of leaf s)
    ! just for reporting
    if (Ed>0.0) then
      w_scale2=min(1.0,ws/Ed)
    else
      w_scale2=1.0
    end if

    ! finally, convert units of stomatal conductance to m/s from mol/(m2 s) by
    ! multiplying it by a volume of a mole of gas
    gs = gs * kR * Tl / p_surf
    !write(899, '(25(E12.4,","))') rad_net,par_net,apot*3600*12,acl*3600*12,Ed

    end associate

  end subroutine gs_leuning


  ! subroutine calc_solarzen(td, latdegrees, cosz, solarelev, solarzen)
  !   ! Calculate solar zenith angle **in radians**
  !   ! From Spitters, C. J. T. (1986), AgForMet 38: 231-242.
  !   implicit none
  !   real, intent(in) :: td             ! day(to minute fraction)
  !   real, intent(in) :: latdegrees     ! latitude in degrees
  !   real :: hour,latrad
  !   real :: delta    ! declination angle
  !   real :: pi, rad
  !   real, intent(out) :: cosz        ! cosz=cos(zen angle)=sin(elev angle)
  !   real, intent(out) :: solarelev    ! solar elevation angle (rad)
  !   real, intent(out) :: solarzen     ! solar zenith angle (rad)

  !   pi  = 3.1415926
  !   rad = pi / 180.0 ! Conversion from degrees to radians.
  !   hour = (td-floor(td))*24.0
  !   latrad = latdegrees*rad
  !   delta  = asin(-sin(rad*23.450)*cos(2.0*pi*(td+10.0)/365.0))
  !   cosz = sin(latrad)*sin(delta) + &
  !   cos(latrad)*cos(delta)*cos(rad* 15.0*(hour-12.0))
  !   cosz = max (cosz, 0.01)  ! Sun's angular is 0.01
  !   ! compute the solar elevation and zenth angles below
  !   solarelev = asin(cosz)/pi*180.0  !since asin(cos(zen))=pi/2-zen=elev
  !   solarzen = 90.0 - solarelev ! pi/2.d0 - solarelev

  ! end subroutine calc_solarzen


  subroutine getpar_modl_gpp()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads module-specific parameters from input file.
    !----------------------------------------------------------------
    ! unit cost of carboxylation
    params_gpp%beta  = 146.000000

    ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
    params_gpp%rd_to_vcmax  = 0.01400000

    ! Apply identical temperature ramp parameter for all PFTs
    params_gpp%tau_acclim     = 30.0
    params_gpp%soilm_thetastar= 0.6 * 250
    params_gpp%soilm_betao    = 0.0

    ! Jmax cost ratio
    params_gpp%kc_jmax  = 0.41

    ! quantum yield efficiency at optimal temperature, phi_0 (Stocker et al., 2020 GMD Eq. 10)
    params_gpp%kphio = 0.05

    ! shape parameter of temperature-dependency of quantum yield efficiency
    params_gpp%kphio_par_a = 0.0

    ! optimal temperature of quantum yield efficiency
    params_gpp%kphio_par_b = 25.0

  end subroutine getpar_modl_gpp

  function qscomp(T, p) result(qsat)
    !--------Output
    real :: qsat ! Output type: saturated specific humidity, kg/kg
    !--------Inputs
    real :: T    ! temperature, degK
    real :: p    ! pressure, Pa
    !--------local var
    real :: myesat ! sat. water vapor pressure
    real :: Temp ! degC

    ! calculate saturated specific humidity
    Temp = T - 273.16 ! degC
    myesat=MIN(610.78*exp(17.27*Temp/(Temp+237.3)), p) ! Pa
    qsat = 0.622*myesat /(p - 0.378*myesat )

  end function qscomp

end module md_gpp_biomee
