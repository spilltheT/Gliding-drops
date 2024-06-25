#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phaseTF.h"
#include "tension.h"
#include "distance.h"
#include "adapt_wavelet_limited_v2.h"

#define MAXlevel 12                                             // maximum level
#define MINlevel 3                                               // maximum level
#define Ldomain 40
#define tmax 60
#define tsnap (0.05)

// Error tolerances
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define VelErr (1e-2) // error tolerances in velocity
#define OmegaErr (1e-3) // error tolerances in vorticity

double Oh1, Oh2, Oh3;
double RHO31;
double BoX; 

#define Xdist (1.040)
#define Ydist (35.0)
#define R2Drop(x,y,z) (sq(x - Xdist) + sq(y - Ydist))

u.t[left] = dirichlet(0.0);
f1[left] = dirichlet(0.0);
f2[left] = dirichlet(1.0);

int main(int argc, char const *argv[]) {
  
  Oh1 = 1.0;
  Oh2 = 0.1;
  Oh3 = 0.00001; // air has a smaller density

  BoX = 1.0;
  
  RHO31 = 0.001;

  L0=Ldomain;
  X0=-0.3; Y0=0.;
  init_grid (1 << (5));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.0; mu1 = Oh1;
  mu2 = Oh2;
  rho3 = RHO31; mu3 = Oh3;

  f1.sigma = 1.0;
  f2.sigma = 0.5;  

  run();
}


event acceleration(i++) {
  face vector av = a;
  foreach_face(x){
    av.x[] -= BoX;
  }
}

event init(t = 0){
  if(!restore (file = "dump")){
    refine((R2Drop(x,y,z) < 1.44) && (level < MAXlevel));
    fraction (f1, 1. - R2Drop(x,y,z));
    fraction (f2, -x);
    foreach () {
      u.x[] = -0.5*f1[];
      u.y[] = -3.0*f1[];
    }
    boundary((scalar *){f1, f2, u.x, u.y});
  }
}

int refRegion(double x, double y, double z){

  return (x > -0.35 && x < 0.1 && y < 37.0 && y > 33.0 ? MAXlevel:
          x > -0.5 && x < 1.80 && y < 40.0 && y > 30.0 ? MAXlevel-1:
          x > -1.0 && x < 2.0 && y < 40.0 && y > 25.0 ? MAXlevel-2:
          x > -2.0 && x < 3.0 && y < 43.0 && y > 20.0 ? MAXlevel-3:
          x < 4.0 && y < 45.0 && y > 10.0 ? MAXlevel-4:
          MAXlevel-5
        );
}

event adapt(i++) {

  scalar KAPPA1[], KAPPA2[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  boundary ((scalar *){KAPPA1, KAPPA2});
  adapt_wavelet_limited ((scalar *){f1, f2, u.x, u.y, KAPPA1, KAPPA2},
    (double[]){fErr, fErr, VelErr, VelErr, KErr, KErr},
    refRegion, MINlevel);
}

event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i+=10) {
  
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += f1[]*(sq(u.x[]) + sq(u.y[]))*(2*pi*y)*sq(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

}
