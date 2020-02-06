# Fortran-code-to-solve-Schrodinger-equation
*Fortran code to solve Schrodinger equation: Solution to Hydrogen & Helium atom by LAPACK subroutines*

The difference among these solvers dsyev, dsbev and dstebz resides in the execution time. In the case of dstebz, the execution time is faster as seen by local CPU_time. Since, in this routine eigenvalues can be
called in the given range, which implies no need for extra calculation. If compared to dsbev subroutine the diagonal and sub or super-diagonal term are stored in another matrix of these dimensions which
is a little faster. Lastly, dsyev just takes the given matrix and solves for eigenvalues and optionally eigenvectors. Over here it doesn’t do any extra effort of assigning matrix to other dimensions and solve.
So, if by execution time (by calling CPU_time subroutine) it could be stated in this way

dstebz < dsbev < dsyev
