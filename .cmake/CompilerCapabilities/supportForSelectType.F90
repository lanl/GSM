
program main
  use, intrinsic :: iso_fortran_env, only: int8, real64

  integer(int8) :: foo
  real(real64) :: bar

  select type (foo)
     type is (int8)
        write(*,"(A)") "'Foo' is an 'int8' variable."
     class default
        error stop "The data type of 'Foo' was NOT recongized correctly."
  end select

  select type (bar)
     type is (real64)
        write(*,"(A)") "'Bar' is a 'real64' variable."
     class default
        error stop "The data type of 'Bar' was NOT recongized correctly."
  end select

end program main
