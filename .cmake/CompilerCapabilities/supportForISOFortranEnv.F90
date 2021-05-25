subroutine foo()
    use, intrinsic:: iso_fortran_env, only: int32, int64, real32, real64, &
        & output_unit, error_unit

    integer(int32) :: stdInt
    integer(int64) :: bigInt
    real(real32)   :: smallReal
    real(real64)   :: stdReal

    write(output_unit, *) "Writes to ouput_unit"
    write(error_unit, *) "Writes to error_unit"
end subroutine foo

program main
end program main
