This is Digital C Library (dcl) created to with the aim to replace old fortran functions and complete rewrite
some of them.

Currently, this library supports the following functions:

1. Gamma function.
I used Lanczos approximation with special parameters oriented to work with double (64 bits) type.
I tested it for various cases and the minimum number of significant digits is 14(for large numbers).
The calculation is very fast.

2. Orthogonal polynoms.
Realised now:
Legendre, Chebyshev, Laguerre and Hermite polynoms. I use the recurrent relations with checks for negative values.
In some functions i use the recurrence directly which may be slower for large orders.

3. Elliptic functions.
EllipticE and EllipticK tested for various values (include very close to 1) and seems to work fast and gives good
accuracy.

4. SphericalHarmonicY
This is complex function based on orthogonal polynoms.

5. ExpIntegralEi
Used series at zero point and asymptotic expansion at infinity. Tested and have good accuracy.

6. GaussIntegrate
Provided simple interface to gauss rule based on LegendrePoly (weights are precomputed). Double and complex.

7. Bessel Functions
Use external fortran realization for complex Bessel functions (Amos library)

8. Linear algebra
QR decomposition and multiplication of matrices.

9. Clebsch Gordan coefficients

10. Incomplete gamma function and ExpIntegral function
I use continued fraction representation which is fast and have good accuracy
