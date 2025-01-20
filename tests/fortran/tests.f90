module tests
    !////////////////////////////////////////////////////////////////
    ! Module containing fortran unit tests
    ! Attention: works in combination with default_insert = .true.
    !----------------------------------------------------------------
    use datatypes_biomee

    implicit none
    private
    public test_linked_list

    contains

    subroutine test_linked_list()

      ! Inout variable
      type(vegn_tile_type) :: vegn
      type(cohort_item), pointer :: ptr
      integer :: start_uid

      if (vegn%n_cohorts() /= 0) STOP 'Wrong cohort number (expected 0)'

      ptr => vegn%new_cohort()
      ptr%cohort%height = 10
      start_uid = ptr%uid - 1

      if (vegn%n_cohorts() /= 1) STOP 'Wrong cohort number (expected 1)'

      ptr => vegn%new_cohort()
      ptr%cohort%height = 5

      if (vegn%n_cohorts() /= 2) STOP 'Wrong cohort number (expected 2)'
      if (ptr%uid - start_uid /= 2) STOP 'Wrong uid (expected: 2)'
      if( vegn%heap%cohort%height /= 5) STOP 'Wrong value (expected 5)'

      ptr => vegn%new_cohort()
      ptr%cohort%height = 6

      if (vegn%n_cohorts() /= 3) STOP 'Wrong cohort number (expected 3)'
      if (ptr%uid - start_uid /= 3) STOP 'Wrong uid (expected: 3)'
      if( vegn%heap%cohort%height /= 6) STOP 'Wrong value (expected 6)'

      ptr => vegn%remove_cohort(2 + start_uid) ! second (midle)

      if (vegn%n_cohorts() /= 2) STOP 'Wrong cohort number (expected 2)'
      if( ptr%uid - start_uid /= 1) STOP 'Wrong uid (expected: 1)'
      if( vegn%heap%cohort%height /= 6) STOP 'Wrong value (expected 6)'
      if( vegn%heap%next%cohort%height /= 10) STOP 'Wrong value (expected 10)'
      if( vegn%heap%next%cohort%height /= 10) STOP 'Wrong value (expected 10)'
      if( ptr%cohort%height /= 10) STOP 'Wrong value (expected 10)'

      if( vegn%heap%uid - start_uid /= 3) STOP 'Wrong uid (expected: 3)'
      ptr => vegn%remove_cohort(3 + start_uid) !first
      if (vegn%n_cohorts() /= 1) STOP 'Wrong cohort number (expected 1)'
      if( vegn%heap%cohort%height /= 10) STOP 'Wrong value (expected 10)'
      if( ptr%cohort%height /= 10) STOP 'Wrong value (expected 10)'

      ptr => vegn%remove_cohort(3 + start_uid) !not existing uid
      if (vegn%n_cohorts() /= 1) STOP 'Wrong cohort number (expected 1)'
      if (associated(ptr)) STOP 'Expected: FALSE'

      ptr => vegn%remove_cohort(1 + start_uid) !last one
      if (vegn%n_cohorts() /= 0) STOP 'Wrong cohort number (expected 0)'
      if (associated(ptr)) STOP 'Expected: FALSE'

      ptr => vegn%new_cohort(.false.)
      ptr%cohort%height = 10

      if (vegn%n_cohorts() /= 1) STOP 'Wrong cohort number (expected 1)'
      if (ptr%uid - start_uid /= 4) STOP 'Wrong uid (expected: 4)'

      ptr => vegn%new_cohort(.false.)
      ptr%cohort%height = 5

      if (vegn%n_cohorts() /= 2) STOP 'Wrong cohort number (expected 2)'
      if (ptr%uid - start_uid /= 5) STOP 'Wrong uid (expected: 5)'

      ptr => vegn%new_cohort(.false.)
      ptr%cohort%height = 2

      if (vegn%n_cohorts() /= 3) STOP 'Wrong cohort number (expected 3)'
      if (ptr%uid - start_uid /= 6) STOP 'Wrong uid (expected: 6)'

      ptr => vegn%remove_cohort(6 + start_uid) !remove tail
      if (vegn%n_cohorts() /= 2) STOP 'Wrong cohort number (expected 2)'
      if (vegn%heap%uid - start_uid /= 4) STOP 'Wrong uid (expected: 4)'
      if (associated(ptr)) STOP 'Expected: FALSE'

      call vegn%clean() !Remove all
      if (vegn%n_cohorts() /= 0) STOP 'Wrong cohort number (expected 0)'

      call vegn%sort_cohorts_by_height(.false.)

      if (vegn%n_cohorts() /= 0) STOP 'Wrong cohort number (expected 0)'

      ptr => vegn%new_cohort()
      ptr%cohort%height = 10

      if (vegn%n_cohorts() /= 1) STOP 'Wrong cohort number (expected 1)'

      call vegn%sort_cohorts_by_height(.false.)

      if (vegn%n_cohorts() /= 1) STOP 'Wrong cohort number (expected 1)'

      ptr => vegn%new_cohort()
      ptr%cohort%height = 2
      ptr => vegn%new_cohort()
      ptr%cohort%height = 5

      if (vegn%n_cohorts() /= 3) STOP 'Wrong cohort number (expected 3)'

      call vegn%sort_cohorts_by_height(.false.)

      if (vegn%n_cohorts() /= 3) STOP 'Wrong cohort number (expected 3)'
      if( vegn%heap%cohort%height /= 10) STOP 'Wrong value (expected 10)'
      if( vegn%heap%next%cohort%height /= 5) STOP 'Wrong value (expected 5)'
      if( vegn%heap%next%next%cohort%height /= 2) STOP 'Wrong value (expected 2)'

      call vegn%sort_cohorts_by_height(.true.)

      if (vegn%n_cohorts() /= 3) STOP 'Wrong cohort number (expected 3)'
      if( vegn%heap%cohort%height /= 2) STOP 'Wrong value (expected 2)'
      if( vegn%heap%next%cohort%height /= 5) STOP 'Wrong value (expected 5)'
      if( vegn%heap%next%next%cohort%height /= 10) STOP 'Wrong value (expected 10)'

      call vegn%clean()

      if (vegn%n_cohorts() /= 0) STOP 'Wrong cohort number (expected 0)'

    end subroutine test_linked_list

end module tests