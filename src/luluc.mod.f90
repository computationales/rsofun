module md_luluc
    !////////////////////////////////////////////////////////////////
    ! Module containing data structures and logic adding Land use &
    ! land use change capabilities to biomee.
    !----------------------------------------------------------------
    implicit none
    private
    public :: update_lu_state, populate_outarray_annual_land_use

    !type luluc_state
    !    real :: fraction
    !end type luluc_state

    contains

    subroutine update_lu_state(lu_state, forcing, n_lu)
    !////////////////////////////////////////////////////////////////
    ! Update the LU fractions 'lu_state' by applying the LU transitions
    ! defined in the square matrix 'forcing'.
    !----------------------------------------------------------------

        ! Arguments
        integer, intent(in) :: n_lu
        real, intent(in), dimension(n_lu, n_lu) :: forcing
        real, intent(inout), dimension(n_lu) :: lu_state

        ! Local variables
        integer :: i, j
        real :: delta

        do i = 1, n_lu
            do j = 1, n_lu
                delta = forcing(i, j)
                lu_state(i) = lu_state(i) - delta
                lu_state(j) = lu_state(j) + delta
            end do
        end do

    end subroutine update_lu_state

    subroutine populate_outarray_annual_land_use(year, lu_state, output_annual_luluc_tile)
        use, intrinsic :: iso_fortran_env, dp=>real64

        ! Arguments
        integer, intent(in)                          :: year
        real, intent(in), dimension(:)               :: lu_state
        real(kind=dp), intent(inout), dimension(:,:) :: output_annual_luluc_tile

        output_annual_luluc_tile(1,:) = dble(year)
        output_annual_luluc_tile(2,:) = dble(lu_state)

    end subroutine populate_outarray_annual_land_use

end module md_luluc