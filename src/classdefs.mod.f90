module md_classdefs
  !////////////////////////////////////////////////////////////////
  ! Module contains Fortran 90 derived-type declarations to define
  ! material pools in SOFUN and functions applicable to pool types.
  ! Pools can be:
  ! - carbon, consisting of C-12
  ! - nitrogen, consisting of N-14
  ! - organic material, consisting of carbon and nitrogen (inherits
  !   their defitions).
  !----------------------------------------------------------------
  use md_params_core, only: eps

  implicit none

  private

  ! Public types
  public carbon, nitrogen, orgpool, common_fluxes

  ! Public functions
  public orgcp, orgcpRec, orgsub, orgmv, &
    orgmvRec, orginit, cmv, cmvRec, ccp, ccpRec, csub, cinit, nmv,  &
    nmvRec, ncp, ncpRec, nsub, ninit, orgfrac, cfrac, nfrac, orgplus, &
    cplus, nplus, orgminus, cminus, nminus, cton, ntoc, update_fluxes


  ! Carbon, so far contains only c12 (to be extended for c13)
  type carbon
   real :: c12 = 0.0
  end type carbon

  type nitrogen
   real :: n14 = 0.0
  end type nitrogen

  ! Organic pools, contain carbon (c12) and nitrogen (n14)
  type orgpool
   type(carbon)   :: c
   type(nitrogen) :: n
  end type orgpool

  type :: common_fluxes
    ! Note: the unit depends on the context
    real    :: Trsp          = 0.0
    real    :: GPP           = 0.0
    real    :: NPP           = 0.0
    real    :: Resp          = 0.0
    real    :: Nup           = 0.0
    real    :: fixedN        = 0.0
  end type common_fluxes

contains
!=========================LOW-LEVEL================================

  subroutine update_fluxes(fluxes, delta)
    ! Add delta quantities to partial fluxes (accounting)
    type(common_fluxes), intent(inout) :: fluxes
    type(common_fluxes), intent(in) :: delta

    fluxes%Trsp   = fluxes%Trsp   + delta%Trsp
    fluxes%GPP    = fluxes%GPP    + delta%GPP
    fluxes%NPP    = fluxes%NPP    + delta%NPP
    fluxes%Resp   = fluxes%Resp   + delta%Resp
    fluxes%Nup    = fluxes%Nup    + delta%Nup
    fluxes%fixedN = fluxes%fixedN + delta%fixedN

  end subroutine update_fluxes

!--------------------------ORGANIC---------------------------------

  subroutine orgcp( amount, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" organic mass to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: to
    real, optional, intent(in) :: scale

    if ( present( scale ) ) then
      call ccp( amount%c,to%c, scale )
      call ncp( amount%n,to%n, scale )
    else
      call ccp( amount%c,to%c )
      call ncp( amount%n,to%n )
    end if

  end subroutine orgcp


  subroutine orgcpRec( amount, to, outc, outn)
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" organic mass to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source.
    !  Additionally records amount copied (adding to outc and outn).
    !----------------------------------------------------------------    
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: to
    real, intent(inout) :: outc
    real, intent(inout) :: outn

    outc = outc + amount%c%c12
    outn = outn + amount%n%n14

    call ccp( amount%c,to%c)
    call ncp( amount%n,to%n)

  end subroutine orgcpRec


  subroutine orgsub( amount, from )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "substract" organic mass ('amount') from source 
    !  pool ('from'). ONLY substracts, but does NOT add.
    !----------------------------------------------------------------
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: from

    call csub( amount%c,from%c)
    call nsub( amount%n,from%n)

  end subroutine orgsub


  subroutine orgmv( amount, from, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" organic mass ('amount') from source pool 
    !  ('from') to destination pool ('to'). Substracts amount moved 
    !  ('amount') from source. 'orgmv' is the combination of 'orgcp' 
    !  and 'orgsub'
    !----------------------------------------------------------------
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: from
    type(orgpool), intent(inout) :: to
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if ( present( scale ) ) then
      call orgcp( orgfrac(scale,amount), to ) 
      call orgsub( amount, from )       
    else
      call orgcp( amount, to )
      call orgsub( amount, from )
    end if  

  end subroutine orgmv


  subroutine orgmvRec( amount, from, to, outc, outn, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" organic mass ('amount') from source pool 
    !  ('from') to destination pool ('to') and to additionally record
    !  amount moved. Substracts amount moved ('amount') from source. 
    ! 'orgmvRec' is the combination of 'orgcp' and 'orgsub'.
    !----------------------------------------------------------------    
    type(orgpool), intent(in) :: amount
    type(orgpool), intent(inout) :: from
    type(orgpool), intent(inout) :: to
    real, intent(inout) :: outc
    real, intent(inout) :: outn
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if ( present( scale ) ) then
      outc = outc + amount%c%c12 * scale
      outn = outn + amount%n%n14 * scale
      call orgcp(orgfrac(scale,amount),to) 
      call orgsub( amount, from )
    else
      outc = outc + amount%c%c12
      outn = outn + amount%n%n14
      call orgcp( amount, to)
      call orgsub( amount, from )
    end if  

  end subroutine orgmvRec


  subroutine orginit( pool )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to initialise organic pool
    !----------------------------------------------------------------
    type(orgpool), intent(inout) :: pool

    call cinit(pool%c)
    call ninit(pool%n)

  end subroutine orginit


  !--------------------------CARBON----------------------------------

  subroutine cmv( amount, from, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" only C from organic mass ('amount') from 
    !  source pool ('from') to destination pool ('to'). Substracts 
    !  amount moved ('amount') from source. 'cmv' is the combination 
    !  of 'ccp' and 'csub'. 
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: from
    type(carbon), intent(inout) :: to
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if ( present( scale ) ) then
      call ccp(cfrac(scale,amount),to)
      call csub( amount, from )
    else
      call ccp( amount, to)
      call csub( amount, from )
    end if

  end subroutine cmv


  subroutine cmvRec( amount, from, to, outc, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" only C from organic mass ('amount') from 
    !  source pool ('from') to destination pool ('to'). Substracts 
    !  amount moved ('amount') from source. 'cmvRec' is the combination 
    !  of 'ccp' and 'csub'. Additionally adds addition to outc/outn.
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: from
    type(carbon), intent(inout) :: to
    real, intent(inout) :: outc
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')


    if ( present( scale ) ) then
      outc = outc + amount%c12 * scale
      call ccp(cfrac(scale,amount),to)
      call csub( amount, from )
    else
      outc = outc + amount%c12
      call ccp( amount, to)
      call csub( amount, from )
    end if

  end subroutine cmvRec


  subroutine ccp( amount, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" carbon to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: to
    real, optional, intent(in) :: scale

    if ( present( scale ) ) then
      to%c12 = to%c12 + amount%c12 * scale
    else
      to%c12 = to%c12 + amount%c12 
    end if

  end subroutine ccp


  subroutine ccpRec( amount, to, outc)
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" carbon to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: to
    real, intent(inout) :: outc

    to%c12 = amount%c12 + to%c12
    outc = outc + amount%c12

  end subroutine ccpRec


  subroutine csub( amount, from )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "substract" organic mass ('amount') from source 
    !  pool ('from'). ONLY substracts, but does NOT add.
    !----------------------------------------------------------------
    type(carbon), intent(in) :: amount
    type(carbon), intent(inout) :: from
    
    from%c12 = from%c12 - amount%c12
     
  end subroutine csub


  subroutine cinit(pool)
    !////////////////////////////////////////////////////////////////
    !  Generic SR to initialise organic pool
    !----------------------------------------------------------------
    type(carbon), intent(inout) :: pool

    pool%c12 = 0.0

  end subroutine cinit


  !--------------------------NITROGEN--------------------------------

  subroutine nmv( amount, from, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" only C from organic mass ('amount') from 
    !  source pool ('from') to destination pool ('to'). Substracts 
    !  amount moved ('amount') from source. 'nmv' is the combination 
    !  of 'ccp' and 'csub'. 
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: from
    type(nitrogen), intent(inout) :: to
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if ( present( scale ) ) then
      call ncp(nfrac(scale,amount),to)
      call nsub( amount, from )
    else
      call ncp( amount, to)
      call nsub( amount, from )
    end if

  end subroutine nmv


  subroutine nmvRec( amount, from, to, outn, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "move" only C from organic mass ('amount') from 
    !  source pool ('from') to destination pool ('to'). Substracts 
    !  amount moved ('amount') from source. 'nmvRec' is the combination 
    !  of 'ccp' and 'csub'. Additionally adds addition to outc/outn.
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: from
    type(nitrogen), intent(inout) :: to
    real, intent(inout) :: outn
    real, intent(in), optional :: scale ! scale source ('from') to be added to destination ('to')

    if ( present( scale ) ) then
      outn = outn + amount%n14 * scale
      call ncp(nfrac(scale,amount),to)
      call nsub( amount, from )
    else
      outn = outn + amount%n14
      call ncp( amount, to)
      call nsub( amount, from )
    end if

  end subroutine nmvRec


  subroutine ncp( amount, to, scale )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" nitrogen to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: to
    real, intent(in), optional :: scale

    if ( present( scale ) ) then
      to%n14 = to%n14 + amount%n14 * scale
    else
      to%n14 = to%n14 + amount%n14 
    end if

  end subroutine ncp


  subroutine ncpRec( amount, to, outn )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "copy" nitrogen to pool (e.g. for output).
    !  Does NOT substract amount moved ('amount') from source
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: to
    real, intent(inout) :: outn

    to%n14 = amount%n14 + to%n14
    outn = outn + amount%n14

  end subroutine ncpRec


  subroutine nsub( amount, from )
    !////////////////////////////////////////////////////////////////
    !  Generic SR to "substract" nitrogen ('amount') from source 
    !  pool ('from'). ONLY substracts, but does NOT add.
    !----------------------------------------------------------------
    type(nitrogen), intent(in) :: amount
    type(nitrogen), intent(inout) :: from
    
    from%n14 = from%n14 - amount%n14

    return

  end subroutine nsub


  subroutine ninit(pool)
    !////////////////////////////////////////////////////////////////
    !  Generic SR to initialise organic pool
    !----------------------------------------------------------------
    type(nitrogen), intent(inout) :: pool

    pool%n14 = 0.0

  end subroutine ninit


  !--------------------------FUNCTIONS--------------------------------

  function orgfrac( frac, from ) result( out_orgfrac )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable of type 'orgpool' and size
    !  of a fraction 'frac' of source pool ('from')
    !----------------------------------------------------------------
    ! arguments
    real, intent(in)          :: frac
    type(orgpool), intent(in) :: from

    ! function return variable
    type(orgpool) :: out_orgfrac

    out_orgfrac%c%c12 = frac * from%c%c12
    out_orgfrac%n%n14 = frac * from%n%n14

  end function orgfrac


  function cfrac( frac, from ) result( out_cfrac )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable of type 'carbon' and size
    !  of a fraction 'frac' of source pool ('from')
    !----------------------------------------------------------------
 
    ! arguments
    real, intent(in) :: frac
    type(carbon), intent(in) :: from

    ! function return variable
    type(carbon) :: out_cfrac

    out_cfrac%c12 = frac * from%c12

  end function cfrac


  function nfrac( frac, from ) result( out_nfrac )
   !////////////////////////////////////////////////////////////////
   !  Generic function to return variable of type 'nitrogen' and size
   !  of a fraction 'frac' of source pool ('from')
   !----------------------------------------------------------------

   ! arguments
   real, intent(in) :: frac
   type(nitrogen), intent(in) :: from
   
   ! function return variable
   type(nitrogen) :: out_nfrac

   out_nfrac%n14 = frac * from%n14

  end function nfrac


  function orgplus( pool1, pool2, pool3, pool4, pool5, pool6, pool7, pool8, pool9, pool10 ) result( out_orgplus )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable sum of two pools of type 
    !  'orgpool'. Sum is of type 'orgpool' as well.
    !----------------------------------------------------------------

    ! arguments
    type(orgpool), intent(in) :: pool1
    type(orgpool), intent(in) :: pool2
    type(orgpool), intent(in),optional :: pool3
    type(orgpool), intent(in),optional :: pool4
    type(orgpool), intent(in),optional :: pool5
    type(orgpool), intent(in),optional :: pool6
    type(orgpool), intent(in),optional :: pool7
    type(orgpool), intent(in),optional :: pool8
    type(orgpool), intent(in),optional :: pool9
    type(orgpool), intent(in),optional :: pool10

    ! function return variable
    type(orgpool) :: out_orgplus

    out_orgplus%c = cplus(pool1%c,pool2%c)
    out_orgplus%n = nplus(pool1%n,pool2%n)

    if (present(pool3)) then
      out_orgplus%c = cplus(out_orgplus%c,pool3%c)
      out_orgplus%n = nplus(out_orgplus%n,pool3%n)
      if (present(pool4)) then
        out_orgplus%c = cplus(out_orgplus%c,pool4%c)
        out_orgplus%n = nplus(out_orgplus%n,pool4%n)
        if (present(pool5)) then
          out_orgplus%c = cplus(out_orgplus%c,pool5%c)
          out_orgplus%n = nplus(out_orgplus%n,pool5%n)
          if (present(pool6)) then
            out_orgplus%c = cplus(out_orgplus%c,pool6%c)
            out_orgplus%n = nplus(out_orgplus%n,pool6%n)
            if (present(pool7)) then
              out_orgplus%c = cplus(out_orgplus%c,pool7%c)
              out_orgplus%n = nplus(out_orgplus%n,pool7%n)
              if (present(pool8)) then
                out_orgplus%c = cplus(out_orgplus%c,pool8%c)
                out_orgplus%n = nplus(out_orgplus%n,pool8%n)
                if (present(pool9)) then
                  out_orgplus%c = cplus(out_orgplus%c,pool9%c)
                  out_orgplus%n = nplus(out_orgplus%n,pool9%n)
                  if (present(pool10)) then
                    out_orgplus%c = cplus(out_orgplus%c,pool10%c)
                    out_orgplus%n = nplus(out_orgplus%n,pool10%n)
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if


  end function orgplus


  function cplus( pool1, pool2, pool3, pool4, pool5, pool6, pool7, pool8, pool9, pool10 ) result( out_cplus )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable sum of two pools of type 
    !  'carbon'. Sum is of type 'carbon' as well.
    !----------------------------------------------------------------

    ! arguments
    type(carbon), intent(in) :: pool1
    type(carbon), intent(in) :: pool2
    type(carbon), intent(in), optional :: pool3
    type(carbon), intent(in), optional :: pool4
    type(carbon), intent(in), optional :: pool5
    type(carbon), intent(in), optional :: pool6
    type(carbon), intent(in), optional :: pool7
    type(carbon), intent(in), optional :: pool8
    type(carbon), intent(in), optional :: pool9
    type(carbon), intent(in), optional :: pool10

    ! function return variable
    type(carbon) :: out_cplus

    out_cplus%c12 = pool1%c12 + pool2%c12

    if (present(pool3)) then
      out_cplus%c12 = out_cplus%c12 + pool3%c12
      if (present(pool4)) then
        out_cplus%c12 = out_cplus%c12 + pool4%c12
        if (present(pool5)) then
          out_cplus%c12 = out_cplus%c12 + pool5%c12
          if (present(pool6)) then
            out_cplus%c12 = out_cplus%c12 + pool6%c12
            if (present(pool7)) then
              out_cplus%c12 = out_cplus%c12 + pool7%c12
              if (present(pool8)) then
                out_cplus%c12 = out_cplus%c12 + pool8%c12
                if (present(pool9)) then
                  out_cplus%c12 = out_cplus%c12 + pool9%c12
                  if (present(pool10)) then
                    out_cplus%c12 = out_cplus%c12 + pool10%c12
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if

  end function cplus


  function nplus( pool1, pool2, pool3, pool4, pool5, pool6, pool7, pool8, pool9, pool10 ) result( out_nplus )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable sum of two pools of type 
    !  'nitrogen'. Sum is of type 'nitrogen' as well.
    !----------------------------------------------------------------

    ! arguments
    type(nitrogen), intent(in) :: pool1
    type(nitrogen), intent(in) :: pool2
    type(nitrogen), intent(in), optional :: pool3
    type(nitrogen), intent(in), optional :: pool4
    type(nitrogen), intent(in), optional :: pool5
    type(nitrogen), intent(in), optional :: pool6
    type(nitrogen), intent(in), optional :: pool7
    type(nitrogen), intent(in), optional :: pool8
    type(nitrogen), intent(in), optional :: pool9
    type(nitrogen), intent(in), optional :: pool10

    ! function return variable
    type(nitrogen) :: out_nplus

    out_nplus%n14 = pool1%n14 + pool2%n14

    if (present(pool3)) then
      out_nplus%n14 = out_nplus%n14 + pool3%n14
      if (present(pool4)) then
        out_nplus%n14 = out_nplus%n14 + pool4%n14
        if (present(pool5)) then
          out_nplus%n14 = out_nplus%n14 + pool5%n14
          if (present(pool6)) then
            out_nplus%n14 = out_nplus%n14 + pool6%n14
            if (present(pool7)) then
              out_nplus%n14 = out_nplus%n14 + pool7%n14
              if (present(pool8)) then
                out_nplus%n14 = out_nplus%n14 + pool8%n14
                if (present(pool9)) then
                  out_nplus%n14 = out_nplus%n14 + pool9%n14
                  if (present(pool10)) then
                    out_nplus%n14 = out_nplus%n14 + pool10%n14
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if

  end function nplus


  function orgminus( pool1, pool2 ) result( out_orgminus )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable difference of two pools of 
    !  type 'orgpool'. Sum is of type 'orgpool' as well.
    !----------------------------------------------------------------

    ! arguments
    type(orgpool), intent(in) :: pool1
    type(orgpool), intent(in) :: pool2

    ! function return variable
    type(orgpool) :: out_orgminus

    out_orgminus%c = cminus( pool1%c, pool2%c )
    out_orgminus%n = nminus( pool1%n, pool2%n )

  end function orgminus


  function cminus( pool1, pool2 ) result( out_cminus )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable difference of two pools of 
    !  type 'carbon'. Sum is of type 'carbon' as well.
    !----------------------------------------------------------------
    ! arguments
    type(carbon), intent(in) :: pool1
    type(carbon), intent(in) :: pool2

    ! function return variable
    type(carbon) :: out_cminus

    out_cminus%c12 = pool1%c12 - pool2%c12

  end function cminus


  function nminus( pool1, pool2 ) result( out_nminus )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return variable difference of two pools of 
    !  type 'carbon'. Sum is of type 'carbon' as well.
    !----------------------------------------------------------------
    ! arguments
    type(nitrogen), intent(in) :: pool1
    type(nitrogen), intent(in) :: pool2

    ! function return variable
    type(nitrogen) :: out_nminus

    out_nminus%n14 = pool1%n14 - pool2%n14

  end function nminus


  function cton( pool, default ) result( out_cton )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return the C:N ratio of an organic pool.
    !----------------------------------------------------------------
    ! arguments
    type(orgpool), intent(in) :: pool
    real, intent(in), optional :: default

    ! function return variable
    real :: out_cton

    if (present(default)) then
      if (pool%n%n14 < eps) then
        out_cton = default
      else
        out_cton = pool%c%c12 / pool%n%n14
      end if
    else
    
    out_cton = pool%c%c12 / pool%n%n14
    
    end if

  end function cton


  function ntoc( pool, default ) result( out_ntoc )
    !////////////////////////////////////////////////////////////////
    !  Generic function to return the N:C ratio of an organic pool.
    !  This is equal to the inverse of the 'cton' function.
    !----------------------------------------------------------------
    ! arguments
    type(orgpool), intent(in) :: pool
    real, intent(in), optional :: default

    ! function return variable
    real :: out_ntoc

    if (present(default)) then
      if (pool%c%c12 < eps) then
        out_ntoc = default
      else
        out_ntoc = pool%n%n14 / pool%c%c12
      end if
    else
    
    out_ntoc = pool%n%n14 / pool%c%c12
    
    end if

  end function ntoc

end module md_classdefs
