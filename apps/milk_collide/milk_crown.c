#include "grid/octree.h" // for 3D simulations
#include "fractions.h"
//#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "SGS.h"  // Eddy viscosity model

#define REDUCED 1
  #if REDUCED
    #include "reduced.h" // reduced gravity = gravitational term applied as interfacial force
  #endif
#include "tension.h"
#include "output_vtu_foreach.h"
// #include "view.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// #include "droplets.h"
// phase 0 is associated with f=0, and rho2, mu2
// phase 1 is associated with f=1 and rho1, mu1
// #define PHASE_OF_INTEREST 1

#define LEVEL 13 // max refinement level
// Reyna's tank parameters
#define L_tank 0.053             // Tank length [m]
#define D_drop 0.0035            // Diameter of BOTH drops [m]
#define xcenter_tank L_tank/2.   // Center of the tank is at
#define zcenter_tank L_tank/2.
#define tank_level 0.0014001     // water level in tank is at y = tank_level
#define GRAVITY 9.81
// harmonic average for viscosity
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

scalar f0[]; // initial volume fraction

//------------------ PARAMS SET IN MAIN ------------------------------
// Droplets positioning
double z_shift, x_shift, x_left, y_left, z_left, x_right, y_right, z_right;
// User input 
double drop_dx, drop_dy, release_height, u0_left, u0_right;
// Parent output directory
char parent_dir[256];
// Numbers characterizing the experiment
double RE_left, RE_right, WE_left, WE_right;


int main(int argc, char* argv[]) {
  /*
    Two drops of equal size shall be released be. Their centers
    are positioned at z = zcenter_tank, x = xcenter_tank \pm xshift and
    y = tank_level + relesase_height + 0 | drop_dy
   */
  
  // Arguments provided are
  // drop_dx = horizontal distance between xcenters of the drops
  // drop_dy = vertical(y) offset right drop center - left drop center
  // release_height = distance of left drop center from milk level
  // Optional
  // u0_left = initial velocity magnitude of the first drop (down)
  // u1_right = initial velocity magnitude of the second drop
  
  if(argc == 6){
    drop_dx = atof(argv[1]); 
    drop_dy = atof(argv[2]);  // This can be negetive
    release_height = atof(argv[3]);
    u0_left = atof(argv[4]); 
    u0_right = atof(argv[5]);
  }
  else{
    assert(argc == 4);

    drop_dx = atof(argv[1]); 
    drop_dy = atof(argv[2]);  // This can be negetive
    release_height = atof(argv[3]);
    
    u0_left = 0.0; 
    u0_right = 0.0;
  }
  // Sanity
  assert(drop_dx > 0); 
  assert(release_height > 0);
  assert(u0_left > -1E-15);
  assert(u0_right > -1E-15);

  // Make sure drops are contained in the domain
  assert(D_drop < L_tank);             // z 
  assert(2*D_drop/2 + drop_dx < L_tank); //[ (   )      (    ) ]
                                         //  <-|<--------->|->
  assert(drop_dx > D_drop); // Don't let them touch initially
  
  assert(release_height+tank_level+D_drop/2 < L_tank);  // 1st is in
  assert(release_height-D_drop/2 > 0);  // Does not touch water

  if(drop_dy > 0){
    assert(release_height+tank_level+D_drop/2 + drop_dy < L_tank);  // So is 2nd
  }
  else{ // Don't touch
    assert(release_height+drop_dy-D_drop/2 > tank_level);
  }

  // Positioning of the drops
  z_shift = 0;
  x_shift = 0.5*(drop_dx - D_drop);
  
  x_left = xcenter_tank - x_shift;
  y_left = release_height+tank_level;
  z_left = zcenter_tank;

  x_right = xcenter_tank + x_shift;
  y_right = y_left + drop_dy;
  z_right = zcenter_tank;

  // physical parameters
  rho1 = 1025.8;
  rho2 =  1.225;
  mu1 = 0.0021;       // water
  mu2 = 1.983E-5;     // air
  const double sigma = 0.0467;
  f.sigma = sigma;
  #if REDUCED
    G.y = -GRAVITY;
  #endif
  // u.x.gradient = u.y.gradient = u.z.gradient = superbee;

  // numerical parameters
  //TOLERANCE = 1e-3;
  //CFL = 0.1;
  DT = 0.1;   // set max time step
  TOLERANCE = 1e-2;
  CFL = 0.4;
  #define femax 0.025 // volume fraction error tolerance
  #define uemax 0.1 // velocity error tolerance

  // grid
  size (L_tank);
  origin (0., 0., 0.);
  init_grid (1 << 5); // might want to initialize at e.g. 1 << (LEVEL-2) due to adaption

  // Set the parent dir
  sprintf(parent_dir,
          "results_dx%g_dy%g_H%g_U0l%g_U0r%g",
          drop_dx, drop_dy, release_height, u0_left, u0_right);
  
  struct stat st = {0};
  if(stat(parent_dir, &st) == -1 && pid() == 0){
      mkdir(parent_dir, 0700);
  }

  // Info on Webber and Reynold numbers based on impact velocity
  const double t_imp_left = (-u0_left + sqrt(pow(u0_left, 2) + 2*GRAVITY*(y_left-tank_level)))/GRAVITY;
  const double U_left = u0_left + t_imp_left*GRAVITY;

  const double t_imp_right = (-u0_right + sqrt(pow(u0_right, 2) + 2*GRAVITY*(y_right-tank_level)))/GRAVITY;
  const double U_right = u0_right + t_imp_right*GRAVITY;

  RE_left = U_left*rho1*D_drop/mu1;
  RE_right = U_right*rho1*D_drop/mu1;

  WE_left = pow(U_left, 2)*rho1*D_drop/sigma;
  WE_right = pow(U_right, 2)*rho1*D_drop/sigma;

  // run
  run();
}


event init (i = 0) {

  // SGS
  Csmag=0.12;

  refine(
   ((fabs(1.-sqrt(pow(x-x_left, 2)+pow(y-y_left, 2)+pow(z-z_left, 2))/(D_drop/2.)) < 0.2) && (level < 10)) ||
   ((fabs(1.-sqrt(pow(x-x_right, 2)+pow(y-y_right, 2)+pow(z-z_right, 2))/(D_drop/2.)) < 0.2) && (level < 10))
  );
  fraction(f0, max(
    max(D_drop/2.-sqrt(pow(x-x_left, 2)+pow(y-y_left, 2)+pow(z-z_left, 2)),
        D_drop/2.-sqrt(pow(x-x_right, 2)+pow(y-y_right, 2)+pow(z-z_right, 2))),
  tank_level-y));

  f0.refine = f0.prolongation = fraction_refine;
  restriction ({f0});

  foreach()
    f[] = f0[];
  boundary({f, f0});

  foreach() {
      u.x[] = 0.;
      u.y[] = 0.;
      u.z[] = 0.;
      // Drop velocities
      if (sqrt(pow(x-x_left, 2) + pow(y-y_left, 2) + pow(z-z_left, 2)) < 1.1*D_drop/2.) {
        u.y[] = -u0_left;
      }

      if (sqrt(pow(x-x_right, 2) + pow(y-y_right, 2) + pow(z-z_right, 2)) < 1.1*D_drop/2.) {
        u.y[] = -u0_right;
      }
  }
}

// print progress to console
event print_progress (i++) {
  if (pid()==0) {
    printf("Iteration %d, t = %g\n", i, t);
    }
}

#if 1
event adapt (i++) {
  //adapt_wavelet ({f}, (double[]){femax}, LEVEL);
  adapt_wavelet ({f, u}, (double[]){femax, uemax, uemax, uemax}, LEVEL);
}
#endif

#if !REDUCED
event acceleration (i++)  {
  face vector av = a;
  foreach_face(y)
    av.y[] = -9.81;
}
#endif

#if 1
event export_vtk (t=0.0001; t += 0.001)
{
  char vtk_dir[512];
  sprintf(vtk_dir, "%s/%s", parent_dir, "vtk");

  struct stat st = {0};
  if (stat(vtk_dir, &st) == -1 && pid() == 0){
    mkdir(vtk_dir, 0700);
  }

  // output the eddy viscosity
  scalar eddy_visc[];
  foreach() {
    eddy_visc[] = 0. ;
    foreach_dimension() {
      eddy_visc[] += (0.5/dimension)*(mu.x[] + mu.x[1]);
    }
  }

  static int j = 0;
  char cmd[1024];
  sprintf (cmd, "%s/snapshot_p%d_%.6i.vtu", vtk_dir, pid(), j++);
  FILE * fp = fopen(cmd, "w");
  //output_vtu_bin_foreach ((scalar *) {f, p, Evis, dt_int, dt_f_int, dt_f_d_int, dt_f_d2_int, dt_f_p_int, dt_f_pr_int}, 
  //                                   (vector *) {u, dt_f_u_int, dt_f_ur_int}, N, fp, false);
  output_vtu_bin_foreach ((scalar *) {f, p}, (vector *) {u}, N, fp, false);
  fclose (fp);
}
#endif

#if 1
// output interface facets as ply file
event output_interface (t=0.0001; t += 0.0005) {
  printf("(left) RE = %g, WE= %g | (right) RE=%g, WE = %g", RE_left, WE_left, RE_right, WE_right);

  char iface_dir[512];
  sprintf(iface_dir, "%s/%s", parent_dir, "interface");

  struct stat st = {0};
  if(stat(iface_dir, &st) == -1){
    mkdir(iface_dir, 0700);
  }

  static int j = 0;
  char name[1024];
  sprintf (name, "%s/interface_p%d_%.6i.vtu", iface_dir, pid(), j++);
  FILE * fp = fopen (name, "w");
  output_vtu (f, fp);

  fclose (fp);
}
#endif

event end (t = 0.002) {
  printf ("Simulation ended");
  printf("(left) RE = %g, WE= %g | (right) RE=%g, WE = %g", RE_left, WE_left, RE_right, WE_right);
}
