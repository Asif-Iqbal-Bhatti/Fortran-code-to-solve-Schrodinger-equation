//-------------------------------- matutil.h --------------------------------
// Contains utility routines for basic operations with vectors and matrices.
// Author: Titus Beu, 2013
//---------------------------------------------------------------------------
#ifndef _MATUTIL_
#define _MATUTIL_

#include <stdio.h>
#include <math.h>
#include "memalloc.h"

//===========================================================================
void MatRead(double **a, int n, int m)
//---------------------------------------------------------------------------
// Reads from the keyboard the elements of matrix a[1:n][1:m]
//---------------------------------------------------------------------------
{
   int i, j;

   for (i=1; i<=n; i++)
      for (j=1; j<=m; j++) { printf("[%i][%i]=",i,j); scanf("%f",&a[i][j]); }
}

//===========================================================================
void MatPrint(double **a, int n, int m)
//---------------------------------------------------------------------------
// Prints the elements of matrix a[1:n][1:m] on the display
//---------------------------------------------------------------------------
{
   int i, j;

   for (i=1; i<=n; i++) {
      for (j=1; j<=m; j++) printf("%12.2e",a[i][j]);
      printf("\n");
   }
}

//===========================================================================
void MatPrintTrans(double **a, int n, int m)
//---------------------------------------------------------------------------
// Prints the elements of transposed matrix a[1:n][1:m] on the display
//---------------------------------------------------------------------------
{
   int i, j;

   for (j=1; j<=m; j++) {
      for (i=1; i<=n; i++) printf("%12.2e",a[i][j]);
      printf("\n");
   }
}

//===========================================================================
void MatZero(double **a, int n, int m)
//---------------------------------------------------------------------------
// Zeros the elements of matrix a[1:n][1:m]
//---------------------------------------------------------------------------
{
   int i, j;

   for (i=1; i<=n; i++)
      for (j=1; j<=m; j++) a[i][j] = 0e0;
}

//===========================================================================
void MatCopy(double **a, double **b, int n, int m)
//---------------------------------------------------------------------------
// Copies matrix a[1:n][1:m] into matrix b[1:n][1:m]
//---------------------------------------------------------------------------
{
   int i, j;

   for (i=1; i<=n; i++) {
      for (j=1; j<=m; j++) b[i][j] = a[i][j];
   }
}

//===========================================================================
void MatTrans(double **a, int n)
//---------------------------------------------------------------------------
// Replaces the square matrix a[1:n][1:n] by its transpose
//---------------------------------------------------------------------------
{
   double t;
   int i, j;

   for (i=2; i<=n; i++)
      for (j=1; j<=(i-1); j++) {
         t = a[i][j]; a[i][j] = a[j][i]; a[j][i] = t;
      }
}

//===========================================================================
void MatDiff(double **a, double **b, double **c, int n, int m)
//---------------------------------------------------------------------------
// Returns the difference of matrices a and b in c[1:n][1:m]
//---------------------------------------------------------------------------
{
   int i, j;

   for (i=1; i<=n; i++)
      for (j=1; j<=m; j++) c[i][j] = a[i][j] - b[i][j];
}

//===========================================================================
void MatProd(double **a, double **b, double **c, int n, int l, int m)
//---------------------------------------------------------------------------
// Returns the product of matrices a[1:n][1:l] and b[1:l][1:m] in c[1:n][1:m]
//---------------------------------------------------------------------------
{
   double t;
   int i, j, k;

   for (i=1; i<=n; i++)
      for (j=1; j<=m; j++) {
         t = 0e0;
         for (k=1; k<=l; k++) t += a[i][k] * b[k][j];
         c[i][j] = t;
      }
}

//===========================================================================
void MatPow(int m, double **a, double **b, int n)
//---------------------------------------------------------------------------
// Returns the m-th power of the square matrix a[1:n][1:n] in b[1:n][1:n]
//---------------------------------------------------------------------------
{
   double **c;
   int k;

   c = Matrix(1,n,1,n);                                         // work array

   MatCopy(a,b,n,n);                                            // case m = 1

   for (k=2; k<=m; k++) {                          // repeated multiplication
      MatProd(a,b,c,n,n,n);
      MatCopy(c,b,n,n);
   }

   FreeMatrix(c,1,1);
}

//===========================================================================
double MatNorm(double **a, int n, int m)
//---------------------------------------------------------------------------
// Returns the max norm of matrix a[1:n][1:m], i.e. max|a[i][j]|
//---------------------------------------------------------------------------
{
   double norm;
   int i, j;

   norm = 0e0;
   for (i=1; i<=n; i++)
      for (j=1; j<=m; j++)
         if (norm < fabs(a[i][j])) norm = fabs(a[i][j]);
   return norm;
}

//===========================================================================
void VecPrint(double a[], int n)
//---------------------------------------------------------------------------
// Prints the elements of vector a[1:n] on the display
//---------------------------------------------------------------------------
{
   int i;

   for (i=1; i<=n; i++) printf("%12.2e",a[i]);
   printf("\n");
}

//===========================================================================
void VecZero(double a[], int n)
//---------------------------------------------------------------------------
// Zeros the elements of vector a[1:n]
//---------------------------------------------------------------------------
{
   int i;

   for (i=1; i<=n; i++) a[i] = 0e0;
}

//===========================================================================
void VecCopy(double a[], double b[], int n)
//---------------------------------------------------------------------------
// Copies vector a[1:n] into vector b[1:n]
//---------------------------------------------------------------------------
{
   int i;

   for (i=1; i<=n; i++) b[i] = a[i];
}

//===========================================================================
void VecDiff(double a[], double b[], double c[], int n)
//---------------------------------------------------------------------------
// Returns the difference of vectors a and b in c[1:n]
//---------------------------------------------------------------------------
{
   int i;

   for (i=1; i<=n; i++) c[i] = a[i] - b[i];
}

//===========================================================================
double VecNorm(double a[], int n)
//---------------------------------------------------------------------------
// Returns the 2-norm of a vector a[1:n]
//---------------------------------------------------------------------------
{
   double norm;
   int i;

   norm = 0e0;
   for (i=1; i<=n; i++) norm += a[i] * a[i];
   norm = sqrt(norm);
   return norm;
}

//===========================================================================
void MatVecProd(double **a, double b[], double c[], int n)
//---------------------------------------------------------------------------
// Returns the product of matrix a[1:n][1:n] and vector b[1:n] in c[1:n]
//---------------------------------------------------------------------------
{
   double t;
   int i, j;

   for (i=1; i<=n; i++) {
      t = 0e0;
      for (j=1; j<=n; j++) t += a[i][j] * b[j];
      c[i] = t;
   }
}

#endif
