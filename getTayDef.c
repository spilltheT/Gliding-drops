#include "navier-stokes/centered.h"
#include "fractions.h"
#include "curvature.h"

char filename[80];

scalar f2[];
double xmax, xmin, ymax, ymin, tayDef;

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  // boundary conditions
  f2[left] = dirichlet(0.0);

  restore (file = filename);
  f2.prolongation = fraction_refine;
  boundary((scalar *){f2, u.x, u.y});

  scalar pos[];
  position (f2, pos, {0,1});

  xmax = statsf(pos).max;
  xmin = statsf(pos).min;

  scalar pos2[];
  position (f2, pos2, {1,0});

  ymax = statsf(pos2).max;
  ymin = statsf(pos2).min;

  tayDef = (xmax - xmin - ymax + ymin)/(xmax - xmin + ymax - ymin);

  boundary((scalar *){f2, u.x, u.y});

  FILE * fp = ferr;
  fprintf (fp, "%g %g %g %g %g\n", t, xmax,xmin, ymax, tayDef);

  fflush (fp);
  fclose (fp);
}
