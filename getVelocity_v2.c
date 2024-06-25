#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80];
double xcm1 , ycm2;

scalar f2[];


int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  // boundary conditions
  f2[left] = dirichlet(0.0);

  restore (file = filename);
  f2.prolongation = fraction_refine;
  boundary((scalar *){f2, u.x, u.y});

  double sumv1 = 0.;
  double sumv2 = 0.;
  double sumf = 0.;

  foreach() {
    sumv1 += clamp(f2[], 0., 1.)*x;
    sumv2 += clamp(f2[], 0., 1.)*y;

    sumf += clamp(f2[], 0., 1.);    
  }
  
//   double xb = 0., yb = 0., vb = 0., sb = 0.;

//   foreach(reduction(+:xb) reduction(+:yb) reduction(+:vb) reduction(+:sb)) {
//     double dv = (1. - f[])*dv();
//     xb += x*dv;
//     yb += y*dv;
//     vb += u.x[]*dv;
//     sb += dv;
//   }
//   area = surface_area(f);

  xcm1 = sumv1/sumf;
  ycm2 = sumv2/sumf;

  boundary((scalar *){f2, u.x, u.y});

  FILE * fp = ferr;
  fprintf(fp, "%f %f %f\n", interpolate(u.x,xcm1, ycm2), interpolate(u.y,xcm1, ycm2), t);

  fflush (fp);
  fclose (fp);
}


// foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
//     double dv = (1. - f[])*dv();
//     xb += x*dv;
//     vb += u.x[]*dv;
//     sb += dv;
//   }
