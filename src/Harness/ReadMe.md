# Harness

The `Harness` library is built to provide a simple interface for which consumers can base a majority of their foundation on by providing a uniform and simple encompassing interface.


## Purpose

The `Harness` library provides several basic features and robust methods by which consumers can rely heavily on. The library is intended to:
+ Provide a central foundation for consumers to use
+ Encompass several robust methods by which consumers can rely on
+ Provide central feature management and incorporation

Some of the features of the `Harness` library include:
+ Parameterized Constants (basic numbers, physical constants, and some emperical data)
+ A lite vector implementation (using a base core object)
+ Lite Fortran exception handling
+ Design-By-Contract protection


## Consuming

The `Harness` library can be consumed by simply importing the module via a Fortran `use` statement. Using the module without the `only` clause will import the entire harness and is the recommended approach for consuming the library.

> The `Harness` library does not attempt to implement all of the features it contains. Instead, it acts as a simple wrapper for the functionality it contains to help simplify the user interface. Details on each individual feature may be found by looking at the associated documentation for that feature. The below code snippet provides an example of the available features:

```fortran
module ConsumerLibrary
    ! Import the Harness library
    !    Note: this example used Constants, Contracts, and a lite Exception
    !          (mostly to demonstrate; Contracts could achieve the same goal)
    use Harness

    implicit none
    
    !> Returns the value of pi (if calc is desired, then randomly calculates based on 100 trials)
    function get_pi(calculate) result(pi_local)
        logical, intent(in), optional :: calculate

        real(real64) :: pi_local = pi

        integer(int16) :: i, hits, numTrials = 100
        real(real64) :: x, y, r

        if (present(calculate) .and. calculate) then
            hits = 0
            do i = 0, numTrials, 1
                ! Randomly calculate pi
                call random_number(x)
                call random_number(y)
                Require(0 <= x .and. x < 1)
                Require(0 <= y .and. y < 1)

                r = ((x**2) + (y**2))^(hlf)
                Validate(zro <= r .and. r <= two, "Calculated R is not within a valid range.")

                if (r < one) then
                    hits += one
                end if
                trails += one
            end do
            pi_local = = four * (real(hits, real64) / real(numTrials, real64))
            if (abs(pi_local - pi) >= 0.2)
                call throw("Warn", "The calculated value of pi is estimated very unprecisely and should not be used for any calcs.")
            end if
        end if
        return
    end function get_pi
end module ConsumerLibrary
```
