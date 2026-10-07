module vilma_precision
   !! Working precision for VILMA2.
   !!
   !! `wp` is double precision and is deliberately kept identical to the C
   !! `double` used by the SHTns and FFTW interfaces, so spectral/spatial arrays
   !! can be passed across the C boundary without copies. The check below makes
   !! that assumption explicit at compile time.
   use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
   use, intrinsic :: iso_c_binding,   only: c_double
   implicit none
   public

   integer, parameter :: sp = real32   !! single precision
   integer, parameter :: dp = real64   !! double precision
   integer, parameter :: wp = real64   !! working precision

   integer, parameter :: i4 = int32
   integer, parameter :: i8 = int64

   ! Fail the build early if the working precision ever stops matching C double,
   ! which would silently break the SHTns/FFTW array passing in vilma_sht. When the
   ! kinds match this is just `real(wp)`; otherwise the kind is -1 (invalid) and
   ! the compiler rejects it.
   real(kind=merge(wp, -1, wp == c_double)), private :: enforce_wp_eq_c_double

end module vilma_precision
