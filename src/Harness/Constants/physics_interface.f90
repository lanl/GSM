!==============================================================================
!> Physical constants interface
!>
!> Provides an interface to arrayed constants for consumer protection
!>
!==============================================================================

    !> Interface to return the value of A^(1/3)
    !>
    !> Interface function for allowing consumers to obtain A^(1/3). This
    !> wrapper provides additional protection to ensure input parameters are not
    !> out of bounds for the internal structure.
    function ato3rd_interface(atomicNumber) result(val)
        !> The atomic number, A, being requested.
        integer(int32), intent(in) :: atomicNumber

        !> The value of A^(1/3).
        real(real64) :: val

        if (atomicNumber < 1) then
            write (error_unit, *) "Cannot have atomic number of 0 or less."
            error stop
        end if
        if (atomicNumber > largest_atomic_weight) then
            write (error_unit, *) "Cannot have atomic number greater than 300."
            error stop
        end if
        val = ato3rd_int(atomicNumber)
        return
    end function ato3rd_interface


