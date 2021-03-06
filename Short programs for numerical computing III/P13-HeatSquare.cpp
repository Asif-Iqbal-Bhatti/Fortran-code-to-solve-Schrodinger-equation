// Steady-state temperature distribution in a thin conducting square plate
#include "memalloc.h"
#include "pde.h"
#include "graphlib.h"
#define pi 3.141592653589793

double Func(double x, double y)                 // RHS function for PoissonXY
   { return 0e0; }
                                // Coefficients for left and right boundaries
void CondX(double y, double &alf_min, double &bet_min, double &gam_min,
                     double &alf_max, double &bet_max, double &gam_max)
{
   alf_min = 1e0; bet_min = 0e0; gam_min = 100e0 * sin(pi*y/10e0);
   alf_max = 1e0; bet_max = 0e0; gam_max = 0e0;
}
                               // Coefficients for lower and upper boundaries
void CondY(double x, double &alf_min, double &bet_min, double &gam_min,
                     double &alf_max, double &bet_max, double &gam_max)
{
   alf_min = 1e0; bet_min = 0e0; gam_min = 0e0;
   alf_max = 1e0; bet_max = 0e0; gam_max = 0e0;
}

int main(int argc, wchar_t** argv)
{
   double **u, *x, *y;
   double eps, hx, hy, umin, umax, xmin, xmax, ymin, ymax;
   int i, j, nx, ny;
   FILE *out;

   xmin = 0e0; xmax = 10e0; ymin = 0e0; ymax = 10e0;     // domain boundaries
   nx = 41; ny = 41;                                 // number of mesh points
   eps = 1e-5;                                 // relative solution tolerance

   u = Matrix(1,nx,1,ny);                                         // solution
   x = Vector(1,nx); y = Vector(1,ny);              // mesh point coordinates

   hx = (xmax-xmin)/(nx-1);
   for (i=1; i<=nx; i++) x[i] = xmin + (i-1)*hx;             // x-mesh points
   hy = (ymax-ymin)/(ny-1);
   for (j=1; j<=ny; j++) y[j] = ymin + (j-1)*hy;             // y-mesh points

   for (j=1; j<=ny; j++)             // initial approximation of the solution
      for (i=1; i<=nx; i++) u[i][j] = 0e0;

   PoissonXY(u,x,y,nx,ny,eps,Func,CondX,CondY);

   out = fopen("Poisson.txt","w");
   fprintf(out,"      x         y          u\n");
   for (j=1; j<=ny; j++)
      for (i=1; i<=nx; i++)
         fprintf(out,"%10.5f%10.5f%15.5e\n",x[i],y[j],u[i][j]);
   fclose(out);

   umin = umax = u[1][1];              // minimum and maximum of the solution
   for (j=1; j<=ny; j++)
      for (i=1; i<=nx; i++) {
         if (u[i][j] < umin) umin = u[i][j];
         if (u[i][j] > umax) umax = u[i][j]; 
      }

   PyGraph w(argc, argv);
   w.GraphInit(800,800);
   w.Contour(u,nx,ny,xmin,xmax,ymin,ymax,umin,umax,
             0.15,0.85,0.15,0.85,"x","y","Plate temperature");
   w.MainLoop();
}
