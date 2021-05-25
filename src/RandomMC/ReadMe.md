# The Random Number Generator Library
___

This random number generator (RNG) was copied from that of MCNP6, code version 2, of the Los Alamos National Laboratory.



## GSM Usage Notes
+ It has been explored to migrate the RNG library to a class-based module. This will NOT work with GSM's use of RNG procedure pointers and the RNG's `rang` function **will** **be** dynamic as a result of being based on the class.
   + If this approach is still desired, **all** of the GSM sub-models **must** include the RNG class and point to their client's object, utilizing the object's `rang` function.
   + The RNG object will need to have embedded within it a procedure pointer to maintain client-control and flexibility, however an `if` statement would need to be included at each instance of obtaining an RNG, reducing the efficiency.
   + This could easily be done, however it would then make each of GSM's sub-models utilizing RNGs dependent on the RNG module