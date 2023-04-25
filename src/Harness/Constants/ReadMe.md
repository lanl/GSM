# Constants

The `Constants` library acts as a simple interface by which consumers can utilize parameterized constants for simulations.


## Purpose

The `Constants` library acts as a method by which consumers can easily use common parameterized numbers. This is generally a very simple task, however abstracting this here provides a few benefits:
+ One-time declaration - consumers do not need to declare their own constants, they can import these with assumed reliability.
+ Robust interface - allows consumers to utilize the module, and no further work is required.
+ Global reach - any adjustments to a constant, like the precision of pi, is propagated throughout all consumers.


## Consuming

The `Constants` library can be easily consumed by those who desire to utilize it. Consumers can import all of the constants, e.g., `use Constants`, or may import only desired constants, e.g., `use Constants, only: <desired constants`.

The below Fortran snippet provides a sample of a consumer whose sole purpose is to return the value of pi. The method below returns the parameterized value by default for runtime efficiency, or if told to calculate pi, will return a very precisely calculated value of pi.

```fortran
module ConsumerLibrary
    ! Import required constants
    use Constants, only: pi

    implicit none
    ...

    !> Calculate Pi, if desired
    function calcPi(calculate) result(pi_local)
        type(logical), intent(in), optional :: calculate

        type(real64) :: pi_local = pi
        if (present(calculate) .and. calculate) then
            pi_local = 4 * atan(1.0)
        end if
        return
    end function calcPi
end module ConsumerLibrary
```
