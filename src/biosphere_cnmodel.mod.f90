module md_biosphere_cnmodel

  use md_params_core
  use md_classdefs
  use md_sofunutils, only: calc_patm
  use md_tile, only: tile_type, tile_fluxes_type, init_tile, init_tile_fluxes, &
    getpar_modl_tile, diag_daily, diag_annual, finalize_tile
  use md_plant, only: getpar_modl_plant
  use md_phenology, only: phenology, getpar_modl_phenology
  use md_waterbal, only: waterbal, solar, getpar_modl_waterbal
  use md_gpp_pmodel, only: getpar_modl_gpp, gpp
  use md_vegdynamics_cnmodel, only: vegdynamics
  use md_soiltemp, only: soiltemp
  use md_npp, only: npp
  use md_nuptake, only: getpar_modl_nuptake, nuptake
  use md_turnover, only: turnover
  use md_landuse, only: landuse
  use md_littersom, only: littersom, getpar_modl_littersom
  use md_ntransform, only: ntransform, getpar_modl_ntransform
  use md_allocation_cnmodel, only: allocation_daily

  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  ! derived types from L1 modules
  type(tile_type),        dimension(nlu) :: tile             ! has gridcell-dimension because values are stored between years
  type(tile_fluxes_type), dimension(nlu) :: tile_fluxes      ! has no gridcell-dimension values need not be recorded

contains

  function biosphere_annual() result( out_biosphere )
    !////////////////////////////////////////////////////////////////
    ! function BIOSPHERE_annual calculates net ecosystem exchange (nee)
    ! in response to environmental boundary conditions (atmospheric 
    ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
    ! LPJ, also formulated as subroutine.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_interface_pmodel, only: myinterface, outtype_biosphere
  
    ! return variable
    type(outtype_biosphere), dimension(ndayyear) :: out_biosphere

    ! local variables
    integer :: dm, moy, doy
    logical, save      :: init_daily            ! is true only on the first day of the simulation 
    logical, parameter :: verbose = .false.     ! change by hand for debugging etc.
    logical, parameter :: baltest = .false.     ! change by hand for debugging etc.
    real               :: cbal1, cbal2, nbal1, nbal2
    type( orgpool )    :: orgtmp1, orgtmp2, orgtmp3, orgtmp4, orgbal1, orgbal2
    real               :: ntmp1, ntmp2, ctmp1, ctmp2

    !----------------------------------------------------------------
    ! INITIALISATIONS
    !----------------------------------------------------------------
    if (myinterface%steering%init) then

      ! set to true on first simulation year and first day
      init_daily = .true.

      !----------------------------------------------------------------
      ! GET MODEL PARAMETERS
      ! read model parameters that may be varied for optimisation
      !----------------------------------------------------------------
      if (verbose) print*, 'getpar_modl() ...'
      call getpar_modl_tile()
      call getpar_modl_plant()
      call getpar_modl_phenology()
      call getpar_modl_waterbal()
      call getpar_modl_gpp()
      call getpar_modl_nuptake()
      call getpar_modl_littersom()
      call getpar_modl_ntransform()
      if (verbose) print*, '... done'

      !----------------------------------------------------------------
      ! Initialise pool variables and/or read from restart file (not implemented)
      !----------------------------------------------------------------
      if (verbose) print*, 'initglobal_() ...'
      call init_tile( tile(:) )
      if (verbose) print*, '... done'

    endif 

    !----------------------------------------------------------------
    ! Get phenology variables (temperature-driven)
    !----------------------------------------------------------------
    if (verbose) print*, 'calling phenology() ...'
    call phenology( tile(:), myinterface%climate(:)%dtemp )
    if (verbose) print*, '... done'

    !----------------------------------------------------------------
    ! LOOP THROUGH MONTHS
    !----------------------------------------------------------------
    doy = 0
    monthloop: do moy=1,nmonth

      !----------------------------------------------------------------
      ! LOOP THROUGH DAYS
      !----------------------------------------------------------------
      dayloop: do dm=1,ndaymonth(moy)
        doy = doy + 1

        if (verbose) print*,'-----------------------------------------------------------'
        if (verbose) print*,'YEAR, DOY ', myinterface%steering%year, doy
        if (verbose) print*,'-----------------------------------------------------------'

        !----------------------------------------------------------------
        ! initialise updated variables (fluxes)
        !----------------------------------------------------------------
        if (verbose) print*,'calling initdaily_() ...'
        call init_tile_fluxes( tile_fluxes(:) )
        if (verbose) print*,'... done.'

        !----------------------------------------------------------------
        ! Get radiation based on daily temperature, sunshine fraction, and 
        ! elevation.
        !----------------------------------------------------------------
        if (verbose) print*,'calling solar() ... '
        if (verbose) print*,'    with argument lat = ', myinterface%grid%lat
        if (verbose) print*,'    with argument elv = ', myinterface%grid%elv
        if (verbose) print*,'    with argument dfsun (ann. mean) = ', sum( myinterface%climate(:)%dfsun / ndayyear )
        if (verbose) print*,'    with argument dppfd (ann. mean) = ', sum( myinterface%climate(:)%dppfd / ndayyear )
        call solar( tile_fluxes(:), &
                    myinterface%grid, & 
                    myinterface%climate(doy),  &
                    doy &
                    )
        if (verbose) print*,'... done'

        !----------------------------------------------------------------
        ! get soil moisture, and runoff
        !----------------------------------------------------------------
        if (verbose) print*,'calling waterbal() ... '
        call waterbal(  tile(:), &
                        tile_fluxes(:), &
                        myinterface%grid, &
                        myinterface%climate(doy) &
                        )
        if (verbose) print*,'... done'

        !----------------------------------------------------------------
        ! calculate soil temperature
        !----------------------------------------------------------------
        if (verbose) print*, 'calling soiltemp() ... '
        call soiltemp( tile(:), &
                       myinterface%climate(:)%dtemp, &
                       moy, &
                       doy, &
                       init_daily &
                       )
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*,'               soiltemp = ', tile(1)%soil%phy%temp
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! update canopy and stand variables and simulate daily 
        ! establishment / sprouting
        !----------------------------------------------------------------
        if (verbose) print*, 'calling vegdynamics() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl 
        if (verbose) print*, '              lai   = ', tile(1)%plant(1)%lai_ind
        !----------------------------------------------------------------
        call vegdynamics( tile(:), tile_fluxes(:), doy, myinterface%steering%init )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl
        if (verbose) print*, '              lai   = ', tile(1)%plant(1)%lai_ind
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! calculate GPP
        !----------------------------------------------------------------
        if (verbose) print*,'calling gpp() ... '
        call gpp( tile(:), &
                  tile_fluxes(:), &
                  myinterface%pco2, &
                  myinterface%climate(doy), &
                  myinterface%climate_memory(doy), &
                  myinterface%vegcover(doy), &
                  myinterface%grid, &
                  myinterface%params_siml%soilmstress, &
                  myinterface%params_siml%tempstress, &
                  init_daily &
                  )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              dgpp  = ', tile_fluxes(1)%plant(1)%dgpp
        if (verbose) print*, '              drd   = ', tile_fluxes(1)%plant(1)%drd
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! substract autotrophic respiration to get NPP, remainder is added 
        ! to labile pool (plabl)
        !----------------------------------------------------------------
        if (verbose) print*, 'calling npp() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              pleaf = ', tile(1)%plant(1)%pleaf
        if (verbose) print*, '              proot = ', tile(1)%plant(1)%proot
        if (verbose) print*, '              psapw = ', tile(1)%plant(1)%psapw
        if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl
        if (verbose) cbal1 = tile(1)%plant(1)%plabl%c%c12
        !----------------------------------------------------------------
        call npp( tile(:), tile_fluxes(:), myinterface%climate(doy) )                  
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              dgpp   = ', tile_fluxes(1)%plant(1)%dgpp
        if (verbose) print*, '              drleaf = ', tile_fluxes(1)%plant(1)%drleaf
        if (verbose) print*, '              drroot = ', tile_fluxes(1)%plant(1)%drroot
        if (verbose) print*, '              drsapw = ', tile_fluxes(1)%plant(1)%drsapw
        if (verbose) print*, '              dnpp   = ', tile_fluxes(1)%plant(1)%dnpp%c12
        if (verbose) print*, '              dcex   = ', tile_fluxes(1)%plant(1)%dcex
        if (verbose) print*, '              plabl  = ', tile(1)%plant(1)%plabl
        if (verbose) print*, '    --- balance: '
        if (verbose) cbal1 = tile_fluxes(1)%plant(1)%dnpp%c12 - tile_fluxes(1)%plant(1)%dcex &
                             - (tile(1)%plant(1)%plabl%c%c12 - cbal1)
        if (verbose) print*, '        gpp - npp - ra_maint          = ', cbal1
        if (baltest .and. abs(cbal1) > eps) stop 'balance 1 not satisfied'
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! calculate N acquisition as a function of C exudation
        !----------------------------------------------------------------
        if (verbose) print*, 'calling nuptake() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              ninorg = ', tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14
        if (verbose) print*, '              nlabl  = ', tile(1)%plant(1)%plabl%n%n14
        if (verbose) ntmp1 = tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14
        if (verbose) ntmp2 = tile(1)%plant(1)%plabl%n%n14
        !----------------------------------------------------------------
        call nuptake( tile(:), tile_fluxes(:) )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              dnup   = ', tile_fluxes(1)%plant(1)%dnup
        if (verbose) print*, '              ninorg = ', tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14
        if (verbose) print*, '              nlabl  = ', tile(1)%plant(1)%plabl%n%n14
        if (verbose) print*, '    --- balance: '
        if (verbose) nbal1 = tile_fluxes(1)%plant(1)%dnup_act + ( tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14 - ntmp1 ) 
        if (verbose) nbal2 = ( tile(1)%plant(1)%plabl%n%n14 - ntmp2 ) + ( tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14 - ntmp1 )
        if (verbose) print*, '        nup - dninorg     = ', nbal1    ! xxx weird result!
        if (verbose) print*, '        dnlabl - dninorg  = ', nbal2
        if (baltest .and. abs(nbal1) > eps) stop 'balance 1 not satisfied'
        if (baltest .and. abs(nbal2) > eps) stop 'balance 2 not satisfied'
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! leaf, sapwood, and fine-root turnover
        !----------------------------------------------------------------
        if (verbose) print*, 'calling turnover() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              lai   = ', tile(1)%plant(1)%lai_ind
        if (verbose) print*, '              pleaf = ', tile(1)%plant(1)%pleaf
        if (verbose) print*, '              proot = ', tile(1)%plant(1)%proot
        if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl
        if (verbose) print*, '              pseed = ', tile(1)%plant(1)%pseed
        if (verbose) print*, '              plitt af = ', tile(1)%soil%plitt_af
        if (verbose) print*, '              plitt as = ', tile(1)%soil%plitt_as
        if (verbose) print*, '              plitt bg = ', tile(1)%soil%plitt_bg
        if (verbose) print*, '              plitt tot = ', orgplus( tile(1)%soil%plitt_af, tile(1)%soil%plitt_as, &
                                                                    tile(1)%soil%plitt_bg )
        if (verbose) orgtmp1 = orgplus( tile(1)%plant(1)%pleaf, &
                                        tile(1)%plant(1)%proot, &
                                        tile(1)%plant(1)%plabl, &
                                        tile(1)%plant(1)%pseed )
        if (verbose) orgtmp2 = orgplus( tile(1)%soil%plitt_af, tile(1)%soil%plitt_as, tile(1)%soil%plitt_bg )
        !----------------------------------------------------------------
        call turnover( tile(:), tile_fluxes(:), doy )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              lai   = ', tile(1)%plant(1)%lai_ind
        if (verbose) print*, '              pleaf = ', tile(1)%plant(1)%pleaf
        if (verbose) print*, '              proot = ', tile(1)%plant(1)%proot
        if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl
        if (verbose) print*, '              pseed = ', tile(1)%plant(1)%pseed
        if (verbose) print*, '              plitt af = ', tile(1)%soil%plitt_af
        if (verbose) print*, '              plitt as = ', tile(1)%soil%plitt_as
        if (verbose) print*, '              plitt bg = ', tile(1)%soil%plitt_bg
        if (verbose) print*, '              plitt = ', orgplus( tile(1)%soil%plitt_af, tile(1)%soil%plitt_as, tile(1)%soil%plitt_bg)
        if (verbose) print*, '   --- balance: '
        if (verbose) orgbal1 = orgminus( orgminus( orgplus( tile(1)%soil%plitt_af, tile(1)%soil%plitt_as, tile(1)%soil%plitt_bg ), &
                                                            orgtmp2   ), orgminus(   orgtmp1,   orgplus( tile(1)%plant(1)%pleaf,   &
                                                            tile(1)%plant(1)%proot, tile(1)%plant(1)%plabl, tile(1)%plant(1)%pseed &
                                                            ) ) )
        if (verbose) print*, '       dlitt - dplant                = ', orgbal1
        if (baltest .and. abs(orgbal1%c%c12) > eps) stop 'balance not satisfied for C'
        if (baltest .and. abs(orgbal1%n%n14) > eps) stop 'balance not satisfied for N'
        if (verbose) print*, '... done'

        ! !----------------------------------------------------------------
        ! ! grass / crop harvest
        ! !----------------------------------------------------------------
        ! if (verbose) print*, 'calling landuse() ... '
        ! if (verbose) print*, '              with state variables:'
        ! if (verbose) print*, '              pleaf = ', tile(1)%plant(1)%pleaf
        ! if (verbose) print*, '              proot = ', tile(1)%plant(1)%proot
        ! if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl
        ! if (verbose) print*, '              mharv = ', tile_fluxes(1)%plant(1)%dharv
        ! if (verbose) orgtmp1 =  orgplus( tile(1)%plant(1)%pleaf, tile(1)%plant(1)%proot, tile(1)%plant(1)%plabl )
        ! if (verbose) orgtmp2 =  tile_fluxes(1)%plant(1)%dharv
        ! !----------------------------------------------------------------
        ! call landuse( tile(:), tile_fluxes(:), doy )
        ! !----------------------------------------------------------------
        ! if (verbose) print*, '              ==> returned: '
        ! if (verbose) print*, '              pleaf = ', tile(1)%plant(1)%pleaf
        ! if (verbose) print*, '              proot = ', tile(1)%plant(1)%proot
        ! if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl
        ! if (verbose) print*, '              mharv = ', tile_fluxes(1)%plant(1)%dharv
        ! if (verbose) print*, '    --- balance: '
        ! if (verbose) orgbal1 = orgminus( orgminus( orgtmp1, orgplus( tile(1)%plant(1)%pleaf, tile(1)%plant(1)%proot, &
        !                         tile(1)%plant(1)%plabl ) ), orgminus( tile_fluxes(1)%plant(1)%dharv, orgtmp2 ) )
        ! if (verbose) print*, '        dharv - dplant  = ', orgbal1
        ! if (baltest .and. abs(orgbal1%c%c12) > eps) stop 'balance not satisfied for C'
        ! if (baltest .and. abs(orgbal1%n%n14) > eps) stop 'balance not satisfied for N'
        ! if (verbose) print*, '... done'

        ! !----------------------------------------------------------------
        ! ! litter and soil decomposition and N mineralisation
        ! !----------------------------------------------------------------
        ! if (verbose) print*, 'calling littersom() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              plitt tot=  ', orgplus( tile(1)%soil%plitt_af, tile(1)%soil%plitt_as, &
                                                                    tile(1)%soil%plitt_bg )
        if (verbose) print*, '              psoil tot = ', orgplus( tile(1)%soil%psoil_fs, tile(1)%soil%psoil_sl )
        if (verbose) print*, '              pexud     = ', tile(1)%soil%pexud
        if (verbose) print*, '              pninorg=    ', tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14
        if (verbose) print*, '              drhet     = ', tile_fluxes(1)%soil%drhet
        if (verbose) print*, '              dnetmin   = ', tile_fluxes(1)%soil%dnetmin
        if (verbose) cbal1 =  tile(1)%soil%plitt_af%c%c12 &
                            + tile(1)%soil%plitt_as%c%c12 &
                            + tile(1)%soil%plitt_bg%c%c12 &
                            + tile(1)%soil%pexud%c12 &
                            + tile(1)%soil%psoil_fs%c%c12 & 
                            + tile(1)%soil%psoil_sl%c%c12 &
                            + tile_fluxes(1)%soil%drhet%c12
        if (verbose) nbal1 =  tile(1)%soil%plitt_af%n%n14 &
                            + tile(1)%soil%plitt_as%n%n14 &
                            + tile(1)%soil%plitt_bg%n%n14 &
                            + tile(1)%soil%psoil_fs%n%n14 & 
                            + tile(1)%soil%psoil_sl%n%n14 &
                            + tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14 &
                            + tile_fluxes(1)%soil%dnfix_free
        !----------------------------------------------------------------
        call littersom( tile(:), tile_fluxes(:), myinterface%climate(doy), doy )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              plitt  = ', orgplus( tile(1)%soil%plitt_af, tile(1)%soil%plitt_as, &
                                                                  tile(1)%soil%plitt_bg )
        if (verbose) print*, '              psoil  = ', orgplus( tile(1)%soil%psoil_fs, tile(1)%soil%psoil_sl )
        if (verbose) print*, '              pninorg= ', tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14
        if (verbose) print*, '              drhet  = ', tile_fluxes(1)%soil%drhet
        if (verbose) print*, '              dnetmin= ', tile_fluxes(1)%soil%dnetmin
        if (verbose) print*, '   --- balance: '
        if (verbose) cbal2 =  tile(1)%soil%plitt_af%c%c12 &
                            + tile(1)%soil%plitt_as%c%c12 &
                            + tile(1)%soil%plitt_bg%c%c12 &
                            + tile(1)%soil%pexud%c12 &
                            + tile(1)%soil%psoil_fs%c%c12 & 
                            + tile(1)%soil%psoil_sl%c%c12 &
                            + tile_fluxes(1)%soil%drhet%c12
        if (verbose) nbal2 =  tile(1)%soil%plitt_af%n%n14 &
                            + tile(1)%soil%plitt_as%n%n14 &
                            + tile(1)%soil%plitt_bg%n%n14 &
                            + tile(1)%soil%psoil_fs%n%n14 & 
                            + tile(1)%soil%psoil_sl%n%n14 &
                            + tile(1)%soil%pnh4%n14 + tile(1)%soil%pno3%n14 &
                            + tile_fluxes(1)%soil%dnfix_free
        if (verbose) cbal1 = cbal2 - cbal1
        if (verbose) nbal1 = nbal2 - nbal1
        if (verbose) print*, '       d( csoil + clitt + cexu + drhet ) = ', cbal1
        ! if (baltest .and. abs(cbal1) > eps .and. (.not. myinterface%steering%do_soilequil)) stop 'balance not satisfied for C'
        if (verbose) print*, '       d( nsoil + nlitt + netmin ) = ', nbal1
        ! if (baltest .and. abs(nbal1) > eps .and. (.not. myinterface%steering%do_soilequil)) stop 'balance not satisfied for N'
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! inorganic soil N dynamics (mass balance test only possible inside module)
        !----------------------------------------------------------------
        if (verbose) print*, 'calling ntransform() ... '
        !----------------------------------------------------------------
        call ntransform( tile(:), tile_fluxes(:), myinterface%landuse(doy), sum(myinterface%climate(:)%dprec), doy )
        !----------------------------------------------------------------
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! allocation of labile pools to biomass
        !----------------------------------------------------------------
        if (verbose) print*, 'calling allocation() ... '
        if (verbose) print*, '              with state variables:'
        if (verbose) print*, '              lai   = ', tile(1)%plant(1)%lai_ind
        if (verbose) print*, '              pleaf = ', tile(1)%plant(1)%pleaf
        if (verbose) print*, '              proot = ', tile(1)%plant(1)%proot 
        if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl
        if (verbose) print*, '              pseed = ', tile(1)%plant(1)%pseed
        if (verbose) print*, '              drgrow= ', tile_fluxes(1)%plant(1)%drgrow
        if (verbose) print*, '              dnup  = ', tile_fluxes(1)%plant(1)%dnup%n14
        if (verbose) cbal1 = tile(1)%plant(1)%pleaf%c%c12 &
                           + tile(1)%plant(1)%proot%c%c12 &
                           + tile(1)%plant(1)%plabl%c%c12 &
                           + tile(1)%plant(1)%pseed%c%c12 &
                           + tile_fluxes(1)%plant(1)%drgrow
        if (verbose) nbal1 = tile(1)%plant(1)%pleaf%n%n14 &
                           + tile(1)%plant(1)%proot%n%n14 &
                           + tile(1)%plant(1)%plabl%n%n14 &
                           + tile(1)%plant(1)%pseed%n%n14 &
                           - tile_fluxes(1)%plant(1)%dnup%n14
        !----------------------------------------------------------------
        call allocation_daily(  tile(:), &
                                tile_fluxes(:),&
                                myinterface%climate(doy), &
                                myinterface%climate_memory(doy) &
                                )
        !----------------------------------------------------------------
        if (verbose) print*, '              ==> returned: '
        if (verbose) print*, '              lai   = ', tile(1)%plant(1)%lai_ind
        if (verbose) print*, '              pleaf = ', tile(1)%plant(1)%pleaf
        if (verbose) print*, '              proot = ', tile(1)%plant(1)%proot
        if (verbose) print*, '              plabl = ', tile(1)%plant(1)%plabl
        if (verbose) print*, '              pseed = ', tile(1)%plant(1)%pseed
        if (verbose) print*, '              drgrow= ', tile_fluxes(1)%plant(1)%drgrow
        if (verbose) print*, '              dnup  = ', tile_fluxes(1)%plant(1)%dnup%n14
        if (verbose) print*, '   --- C balance: '
        if (verbose) cbal2 = tile(1)%plant(1)%pleaf%c%c12 &
                           + tile(1)%plant(1)%proot%c%c12 &
                           + tile(1)%plant(1)%plabl%c%c12 &
                           + tile(1)%plant(1)%pseed%c%c12 &
                           + tile_fluxes(1)%plant(1)%drgrow
        if (verbose) nbal2 = tile(1)%plant(1)%pleaf%n%n14 &
                           + tile(1)%plant(1)%proot%n%n14 &
                           + tile(1)%plant(1)%plabl%n%n14 &
                           + tile(1)%plant(1)%pseed%n%n14 &
                           - tile_fluxes(1)%plant(1)%dnup%n14
        if (verbose) cbal1 = cbal2 - cbal1
        if (verbose) nbal1 = nbal2 - nbal1
        if (verbose) print*, '       d( cleaf + croot + clabl + cseed + rgrowth ) =', cbal1
        ! if (baltest .and. abs(cbal1) > eps) stop 'balance not satisfied for C'
        if (verbose) print*, '       d( nleaf + nroot + nlabl + nseed - nup ) =', nbal1
        ! if (baltest .and. abs(nbal1) > eps) stop 'balance not satisfied for C'
        if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! daily diagnostics (e.g., sum over plant within canopy)
        !----------------------------------------------------------------
        call diag_daily( tile(:), tile_fluxes(:), out_biosphere(doy) )

        init_daily = .false.

      end do dayloop

    end do monthloop

    !----------------------------------------------------------------
    ! annual diagnostics
    !----------------------------------------------------------------
    call diag_annual( tile(:), tile_fluxes(:) )

    !----------------------------------------------------------------
    ! close (experimental) files
    !----------------------------------------------------------------
    if (myinterface%steering%finalize) then
      call finalize_tile()
    end if
    
    if (verbose) print*, 'Done with biosphere for this year. Guete Rutsch!'

  end function biosphere_annual

end module md_biosphere_cnmodel