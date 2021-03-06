# Angular motion of a nonlinear pendulum by the Runge-Kutta method
#    u" = -g/l * sin(u) - k * u',   u(0) = u0, u'(0) = u0'
from math import *
from ode import *
from integral import *
from graphlib import *

g = 9.81e0                                       # gravitational acceleration

def Func(t, u, f):                                    # RHS of 1st order ODEs
   f[1] =  u[2]                                         # u[1] = u, u[2] = u'
   f[2] = -g/l * sin(u[1]) - k * u[2]

#============================================================================
def Kel(m):
#----------------------------------------------------------------------------
#  Returns the complete elliptic integral of the 1st kind
#  Calls: qImprop2 (integral.py)
#----------------------------------------------------------------------------
   eps = 1e-7                                            # relative precision
   def fKel(z): return 1e0 / sqrt((1e0-z*z)*(1e0-m*z*z))          # integrand
   return qImprop2(fKel,0e0,1e0,eps)

# main

nn  = [0]*3                                         # ending indexes of plots
col = [""]*3                                                # colors of plots
sty = [0]*3                                                 # styles of plots

l = 1e0                                                     # pendulum length
k = 0e0                                                # velocity coefficient
du0 = 0e0                                                # initial derivative
tmax = 20e0                                                       # time span
ht = 0.001e0                                                 # time step size
u0min = 0.1e0                                                    # minimum u0
u0max = 3.1e0                                                    # maximum u0
hu = 0.1e0                                                 # increment for u0

n = 2                                              # number of 1st order ODEs
nt = int(tmax/ht + 0.5) + 1                            # number of time steps
nu = int((u0max-u0min)/hu + 0.5) + 1                          # number of u0s
u = [0]*(n+1)                                           # solution components
tt = [0]*(2*nt+1); ut = [0]*(2*nt+1)                # arrays for t-dep. plots
tu = [0]*(2*nu+1); uu = [0]*(2*nu+1)               # arrays for u0-dep. plots

for iu in range(1,nu+1):
   u0 = u0min + (iu-1)*hu                              # initial displacement

   t = 0e0; it = 1
   u[1] = u0; u[2] = du0                                     # initial values
   if (iu == 1 ): tt[1   ] = t; ut[1   ] = u[1]             # for smallest u0
   if (iu == nu): tt[1+nt] = t; ut[1+nt] = u[1]              # for largest u0

   nT = 0                                            # number of half-periods
   t1 = t2 = 0e0                                    # bounding solution zeros
   us = u[1]                                                  # save solution
   while (t+ht <= tmax):                                   # propagation loop
      RungeKutta(t,ht,u,n,Func)
      t += ht; it += 1

      if (u[1]*us < 0e0):              # count solution passages through zero
         if (t1 == 0): t1 = t                                  # initial zero
         else: t2 = t; nT += 1                                   # final zero
      us = u[1]                                               # save solution
                                                         # store for plotting
      if (iu == 1 ): tt[it   ] = t; ut[it   ] = u[1]        # for smallest u0
      if (iu == nu): tt[it+nt] = t; ut[it+nt] = u[1]         # for largest u0

   T = 2e0*(t2-t1) / nT                                   # calculated period
   T0 = 2e0*pi*sqrt(l/g)                                    # harmonic period
   Tex = 2/pi * T0 * Kel(pow(sin(0.5e0*u0),2))                 # exact period

   uu[iu   ] = u0; tu[iu   ] = T/T0
   uu[iu+nu] = u0; tu[iu+nu] = Tex/T0

GraphInit(1200,600)

nn[1] = it   ; col[1] = "blue"; sty[1] = 1                  # for smallest u0
nn[2] = it+nt; col[2] = "red" ; sty[2] = 1                   # for largest u0
MultiPlot(tt,ut,ut,nn,col,sty,2,10,0e0,0e0,0,0e0,0e0,0,0.10,0.45,0.15,0.85, \
          "t (s)","u (rad)","Displacement of nonlinear pendulum")

nn[1] =   nu; col[1] = "red" ; sty[1] = 0                # calculated periods
nn[2] = 2*nu; col[2] = "blue"; sty[2] = 1                     # exact periods
MultiPlot(uu,tu,tu,nn,col,sty,2,10,0e0,0e0,0,0e0,0e0,0,0.60,0.95,0.15,0.85, \
          "u0 (rad)","T/T0","Period of nonlinear pendulum")

MainLoop()
