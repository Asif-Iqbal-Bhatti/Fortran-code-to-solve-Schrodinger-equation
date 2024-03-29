# Solves matrix equation by the Gauss-Jordan method
from linsys import *
from matutil import *

n = 4                                                       # order of system
m = 1                                            # number of constant vectors
a = [[0]*(n+1) for _ in range(n+1)]
b = [[0]*(m+1) for _ in range(n+1)]
c = [[0]*(n+1) for _ in range(n+1)]
d = [[0]*(n+1) for _ in range(n+1)]

a[1][1] = 1
a[1][2] = 2
a[1][3] = 3
a[1][4] = 4
b[1][1] = 30
a[2][1] = 2
a[2][2] = 1
a[2][3] = 2
a[2][4] = 3
b[2][1] = 22
a[3][1] = 3
a[3][2] = 2
a[3][3] = 1
a[3][4] = 2
b[3][1] = 18
a[4][1] = 4
a[4][2] = 3
a[4][3] = 2
a[4][4] = 1
b[4][1] = 20
print("A:")
MatPrint(a,n,n)
print("B:")
MatPrint(b,n,m)

MatCopy(a,c,n,n)                                         # save system matrix

det = GaussJordan(a,b,n,m)                       # solve system, inverse in a

print("det A = ",det)
print("Solution:")
MatPrint(b,n,m)

print("Check A^(-1)A = I:")
MatProd(a,c,d,n,n,n)                         # multiply inverse with original
MatPrint(d,n,n)
