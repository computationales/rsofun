module md_cohort_linked_list
  !////////////////////////////////////////////////////////////////
  ! Module containing cohort linked list
  !----------------------------------------------------------------
  use md_cohort

  ! define data types and constants
  implicit none
  private
  !=============== Public types ===========================================================
  public :: cohort_item

  !=============== Public procedures ===========================================================
  public :: next_uid

  integer, private :: CurrentCohortUid = 0

  type :: cohort_item
    integer :: uid ! Unique id. To be filled with 'new_uid()'
    type(cohort_type) :: cohort
    type(cohort_item), pointer :: next => NULL() ! Pointer to next cohort_item. Important to nullify here!

    contains

    procedure has_next

  end type cohort_item

contains

  function has_next(self) result(res)
    logical :: res
    class(cohort_item) :: self

    res = associated(self%next)
  end function has_next

  function next_uid() result(res)
  integer :: res

  CurrentCohortUid = CurrentCohortUid + 1
  res = CurrentCohortUid
end function next_uid

end module md_cohort_linked_list
