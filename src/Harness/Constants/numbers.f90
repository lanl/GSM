!==============================================================================
!> Simple numbers
!>
!> Provides a list of simple numbers, e.g., 1, 2, 1/3, pi, etc., for consumers.
!>
!==============================================================================
  ! Frequently used real numbers:
  real(real64), public, parameter ::  &
       &  zro              = 0.0_real64,                 &
       &  one              = 1.0_real64,                 &
       &  two              = 2.0_real64,                 &
       &  thr              = 3.0_real64,                 &
       &  four             = 4.0_real64,                 &
       &  fiv              = 5.0_real64,                 &
       &  six              = 6.0_real64,                 &
       &  seven            = 7.0_real64,                 &
       &  eight            = 8.0_real64,                 &
       &  nine             = 9.0_real64,                 &
       &  ten              = 10.0_real64,                &
       &  twelve           = 12.0_real64,                &
       &  fourteen         = 14.0_real64,                &
       &  fifteen          = 15.0_real64,                &
       &  twenty           = 20.0_real64,                &
       &  thirty           = 30.0_real64,                &
       &  forty            = 40.0_real64,                &
       &  fortyfiv         = 45.0_real64,                &
       &  fifty            = 50.0_real64,                &
       &  sixty            = 60.0_real64,                &
       &  sxtysix          = 66.0_real64,                &
       &  seventy          = 70.0_real64,                &
       &  eighty           = 80.0_real64,                &
       &  ninety           = 90.0_real64,                &
       &  hundred          = 100.0_real64,               &
       &  one_fifty        = 150.0_real64,               &
       &  one_eighty       = 180.0_real64,               &
       &  two_hundred      = 200.0_real64,               &
       &  three_sixty      = 360.0_real64,               &
       &  thousand         = 1000.0_real64,              &
       &  thsn             = 1000.0_real64,              &
       &  hundred_thousand = 100000.0_real64,            &
       &  million          = 1000000.0_real64,           &
       &  thrd             = one / thr,                  &
       &  hlf              = 0.5_real64,                 &
       &  eighth           = one / eight,                &
       &  twthrd           = two / thr,                  &
       &  tenth            = 0.1_real64,                 &
       &  hundredth        = 0.01_real64,                &
       &  thousandth       = 0.001_real64,               &
       &  millionth        = 1.0e-06_real64,             &
       &  trillionth       = 1.0e-12_real64,             &
       &  pi               = 3.141592653589793d0,        &
                     != pi.  3.14159265358979324_real64,
       &  twpi             = two * pi,                   &
       &  natural_num      = exp(one),                   &
       &  euler            = .577215664901532861_real64, & != Euler constant.
       &  LnTwo            = .693147180559945309_real64, & != Natural logarithm of 2.
       &  SqrtTwo          = 1.41421356237309505_real64    != Square root of 2.


