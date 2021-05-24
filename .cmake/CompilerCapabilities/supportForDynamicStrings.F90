subroutine foo(foodIn)
   character(*), intent(in) :: foodIn
   character(:), allocatable :: string
   character(*), parameter :: pizza = "Pizza"

   string = "Steak"
   string = "Burgers and Fries"
   string = pizza
   write(*, "(A)") string

   write(string, "(A, f8.5, A)") "Pi is ", 3.14159, ", but it's also tasty."
   write(*, "(A)") foodIn, " is a good one too!"

end subroutine foo

program main
    character(:), allocatable:: food
    food = "orange"
    call foo(food)
end program main
