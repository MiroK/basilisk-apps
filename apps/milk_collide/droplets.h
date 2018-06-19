#include "tag.h" // to label individual droplets

// set the phase of interest (0 or 1)


#if PHASE_OF_INTEREST
  #define fi f[]
#else
  #define fi (1.-f[])
#endif


// a struct for symmetric tensors
typedef struct {
  double xx, xy, xz, yy, yz, zz;
} symtens;


scalar dt_int[]; // the integral of dt in each cell (the time!)
scalar dt_f_int[]; // the integral of dt times volume fraction in each cell (the time occupied by f=1)
scalar dt_f_d_int[]; // the integral of dt times diameter (phase 1)
scalar dt_f_d2_int[]; // the integral of dt times squared diameter (phase 1)
vector dt_f_u_int[]; // the integral of dt times velocity(phase 1)
vector dt_f_ur_int[]; // the integral of dt times velocity (phase 2)
scalar dt_f_p_int[]; // the integral of dt times pressure (phase 1)
scalar dt_f_pr_int[]; // the integral of dt times pressure (phase 2)


event initialize_droplets (i=0)
  {

  foreach(){
    dt_int[] = 0.;
    dt_f_int[] = 0.;
    dt_f_d_int[] = 0.;
    dt_f_d2_int[] = 0.;
    dt_f_u_int.x[] = dt_f_u_int.y[] = dt_f_u_int.z[] = 0.;
    dt_f_ur_int.x[] = dt_f_ur_int.y[] = dt_f_ur_int.z[] = 0.;
    dt_f_p_int[] = 0.;
    dt_f_pr_int[] = 0.;
    }
  }


event droplets (i++)
{


  foreach(){
    dt_int[] += dt;
    dt_f_int[] += dt*fi;
    dt_f_p_int[] += dt*fi*p[];
    dt_f_pr_int[] += dt*(1.-fi)*p[];
    foreach_dimension(){
      dt_f_u_int.x[] += dt*fi*u.x[];
      dt_f_ur_int.x[] += dt*(1.-fi)*u.x[];
    }
  }



  
// extract individual droplet information
  scalar m[];
  foreach()
    m[] = fi > 1e-3;
  int n = tag (m);

  /**
  Once each cell is tagged with a unique droplet index, we can easily
  compute the volume *v* and position *b* of each droplet. Note that
  we use *foreach_leaf()* rather than *foreach()* to avoid doing a
  parallel traversal when using OpenMP. This is because we don't have
  reduction operations for the *v* and *b* arrays (yet). */

  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++) // loop over all droplets
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*fi; // volume
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*fi*p.x; // volume integrated position of droplet n is in b[n]
    }

 /**
 When using MPI we need to perform a global reduction to get the
 volumes and positions of droplets which span multiple processes. */

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

// Do another loop and compute shape tensor, radius and squared radius

  symtens shape_tensor[n];

  for (int j = 0; j < n; j++) // loop over all droplets
    shape_tensor[j].xx = shape_tensor[j].xy = shape_tensor[j].xz = shape_tensor[j].yy = shape_tensor[j].yz = shape_tensor[j].zz = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      dt_f_d_int[] += dt*fi*pow(2.54647908947*v[j], 1./3.); // (V*8/pi)**(1/3)
      dt_f_d2_int[] += dt*fi*pow(2.54647908947*v[j], 2./3.); // (V*8/pi)**(2/3), to compute the standard deviation of droplet diameter
      shape_tensor[j].xx += dv()*fi*(x-b[j].x)*(x-b[j].x);
      shape_tensor[j].xy += dv()*fi*(x-b[j].x)*(y-b[j].y);
      shape_tensor[j].xz += dv()*fi*(x-b[j].x)*(x-b[j].z);
      shape_tensor[j].yy += dv()*fi*(x-b[j].y)*(x-b[j].y);
      shape_tensor[j].yz += dv()*fi*(x-b[j].y)*(x-b[j].z);
      shape_tensor[j].zz += dv()*fi*(x-b[j].z)*(x-b[j].z);
    }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, shape_tensor, 6*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif




#if 1
  /**
  Finally we output the volume, position and shape tensor of each droplet to
  standard output. */

  for (int j = 0; j < n; j++)
    // iteration, time, droplet number, volume, position
    fprintf (fout, "%d %g %d %g %g %g %g %g %g %g %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j], shape_tensor[j].xx, shape_tensor[j].xy, shape_tensor[j].xz, shape_tensor[j].yy, shape_tensor[j].yz, shape_tensor[j].zz);
  fflush (fout);
#endif

}

