/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Frank Uhlig (ICP Stuttgart)
   Sources: Jiri Kolafa "Time-Reversible Always Stable Predictor–Corrector
            Method for Molecular Dynamics of Polarizable Molecules"
            J Comput Chem 25: 335-342, 2004
            Infrastructure for memory allocation/deallocation/communicating
            taken from fix_minimize
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_aspc.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "compute.h"
#include "neighbor.h"
#include "domain.h"
#include "modify.h"
#include "pair.h"
#include "output.h"
#include "thermo.h"
#include "timer.h"
#include "error.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixASPC::FixASPC(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg)
{

  if (narg < 5) error->all(FLERR,"Illegal fix ASPC command");

  length = force->inumeric(FLERR,arg[3]);
  nord = length + 2;
  damp = force->numeric(FLERR,arg[4]);

  if ( damp < 0 )
      damp = (length + 2.) / (2.*length + 3.);

  /*FUDO| set some default for what and dim */
  what = NONE;
  dim = 0;

  //FU| some defaults that can be changed by fix_modify
  //FUDO| check and see if we should actually still use the dmax

  //FUDO| actually use dmax?!
  // dmax = 0.01;

  eflag = 1;
  vflag = 0;

  //FUDO| save this to determine the offset later one (we might/should not be inserting here in the first step)
  //FUDO| is this correct here, or do we need to do this in init?
  first = update->ntimestep;

  //FUDO| this is probably dangerous to set to 1 by default
  zerovels = 1;

  nvector = 0;
  peratom = NULL;
  vectors = NULL;
  coeffs = NULL;

  //FUDO| we'll skip the first step, because we have a history of all zeros
  //FUDO| alternatively we could fill up the history with the current displacement
  nonlyhist = 0;
  filledhist = 0;

  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixASPC::~FixASPC()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored data

  memory->destroy(peratom);
  for (int m = 0; m < nvector; m++) memory->destroy(vectors[m]);
  memory->sfree(vectors);
  memory->sfree(coeffs);
}

/* ---------------------------------------------------------------------- */

void FixASPC::init()
{
   //FUDO| set external_force_clear correctly
   external_force_clear = 0;
   
   torqueflag = extraflag = 0;
   if (atom->torque_flag) torqueflag = 1;
   if (atom->avec->forceclearflag) extraflag = 1;

   // allow pair and Kspace compute() to be turned off via modify flags
   // FUDO| this should be included in some way, not sure how
   
   if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
   else pair_compute_flag = 0;
   if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
   else kspace_compute_flag = 0;

   // FUDO| this is prolly totally irrelevant for us
   // orthogonal vs triclinic simulation box
   
   triclinic = domain->triclinic;
   
   bigint ndofme = 3 * static_cast<bigint>(atom->nlocal);
   MPI_Allreduce(&ndofme,&ndoftotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

   // FUDO| perform add_vector here for order vectors with nper dimension
   // FUDO| instead of add_vector just use the general memor allocation routines, because we know what we need

   if ( what == NONE )
      error->all(FLERR, "Don't know what to predict/correct");

   for (int k=0; k<nord; k++)
     add_vector(dim);

   tstart = update->ntimestep;

   generate_coefficients(length);

   // FUDO| this would be nice to do, but we cannot, because we're not in reduced coordinates yet
   // if (nonlyhist == 1)
   //   for ( int n=0; n<nord; n++ )
   //     update_history();

}

void FixASPC::force_clear()
{
  //FUDO| remove?!?
  if (external_force_clear) return;

  // clear global force array
  // if either newton flag is set, also include ghosts

  size_t nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;

  if (nbytes) {
    memset(&atom->f[0][0],0,3*nbytes);
    if (torqueflag) memset(&atom->torque[0][0],0,3*nbytes);
    if (extraflag) atom->avec->force_clear(0,nbytes);
  }
}

void FixASPC::reset_vectors()
{

    // FUDO| here we'll have a second array that will be changed each time step to point to the correct per-atom vectors in the correct length

    nvec = dim * atom->nlocal;

    if ( nvec ) {
        switch ( what ) {
            case COORDS:
                qty = atom->x[0];
                break;
            case QEQ:
                qty = atom->q;
                break;
            case QEQREAX:
                qty = atom->q;
                break;
            default:
                error->all(FLERR, "Do not know which quantity to predict/correct");
                break;
        }
    }

}

void FixASPC::update_history()
{
    //FUDO | copy new history element into vectors

    int i, k;
    double *old, *ndt;

    for ( k=nord-1; k>0; k-- ) {

      old = request_vector(k);
      ndt = request_vector(k-1);

      for ( i=0; i<nvec; i++ )
        old[i] = ndt[i];
    }

    cpy2hist();

}

void FixASPC::cpy2hist()
{
    int i, k;
    double *data;

    data = request_vector(0);

    for ( i=0; i<nvec; i++ )
      data[i] = qty[i];

}

int FixASPC::modify_param(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal fix_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dmax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      dmax = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"zerovels") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1], "yes") == 0) zerovels = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) zerovels = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"kspace_compute_flag") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1], "yes") == 0) kspace_compute_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) kspace_compute_flag = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"pair_compute_flag") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1], "yes") == 0) pair_compute_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) pair_compute_flag = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"damp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      damp = force->numeric(FLERR, arg[iarg+1]);
      if ( damp < 0. )
         damp = (length + 2.) / (2.*length + 3.);
      iarg += 2;
    } else if (strcmp(arg[iarg],"nonlyhist") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      nonlyhist = force->inumeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

  return 2;
}

//FUDO| check in general what to do with the velocties and where to do it
//FUDO| i.e., if we to predict/correct coordinates it would be useful not to additionally integrate it in the verlet cycle, hence zero-ing out the velocities makes sense
//FUDO| this could be done in pre-/post-force and similar

void FixASPC::final_integrate()
{
    reset_vectors();
    update_history();
}

void FixASPC::initial_integrate(int vflag)
{

    //FUDO| need to set the sorted vectors
    //FUOD| not sure, but reset should prolly not be needed

    reset_vectors();

    double **v = atom->v;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;

    if ( zerovels )
      for ( int i=0; i<nlocal; i++ )
        if ( mask[i] & groupbit )
          v[i][0] = v[i][1] = v[i][2] = 0.0;

    if ( !(filledhist) ) {
      // for ( int n=0; n<nord; n++ )
      update_history();

      filledhist = 1;
    }

    if ( ( update->ntimestep - tstart ) < (nonlyhist))
        return;

    predict();
}

void FixASPC::setup_pre_force(int vflag)
{
    //FUDO| there's nothing in the history yet, and nothing to predict from, hence only to the correction
    pre_force(vflag);
    comm->forward_comm();

    update_history();
}

void FixASPC::pre_force(int vflag)
{
    if ( ( update->ntimestep - tstart ) < (nonlyhist))
        return;

    reset_vectors();
    correct();

    //FUDO| we'll not do re-neighboring, should we? 
    comm->forward_comm();
}

void FixASPC::reset_history()
{
    int i, k, baseind;
    int *mask = atom->mask;

}

void FixASPC::predict()
{

    int i, k, n, baseind;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;
    double *data;

    data = request_vector(0);
    baseind = 0;

    int my_nord;
    int tlength = update->ntimestep - tstart;

    if ( tlength <= nord ) {
        my_nord = tlength;
        generate_coefficients(my_nord - 2);
        // printf("TEMPORARY: %i\n", tlength);
        // for ( n=0; n < my_nord; n++ )
        //     printf("COEFF: %14.8f\n", coeffs[n]);
    }
    else
        my_nord = nord;

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit )
        for ( k=0; k<dim; k++) {
           // printf("OLD: %14.8f %14.8f %14.8f\n", data[baseind], data[baseind+1], data[baseind+2]);
           qty[baseind+k] = coeffs[0] * data[baseind+k];
        }

      baseind += dim;
    }

    for ( n=1; n<my_nord; n++ ) {

      baseind = 0;
      data = request_vector(n);

      for ( i=0; i<nlocal; i++ ) {
        if ( mask[i] & groupbit )
          for ( k=0; k<dim; k++)
             qty[baseind+k] += coeffs[n] * data[baseind+k];

        baseind += dim;
      }
    }

}

//FUDO| specific routine needs to be added for each fix
//FUDO| we need a function to provide the scf step that we can then use to perform the damped propagation with
void FixASPC::correct()
{
}

void FixASPC::clear_non_group()
{
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double **f = atom->f;

    for ( int i=0; i<nlocal; i++ )
      if ( !( mask[i] & groupbit ) )
          f[i][0] = f[i][1] = f[i][2] = 0.;
}

int FixASPC::setmask()
{
  int mask = 0;
  // FUDO| what is the correct setup function for INITIAL_INTEGRATE
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

void FixASPC::generate_coefficients(int len)
{
    //FUDO| issue with long long to double casting? do we need that explicitly? 
    // B_j = (-1)^(j+1) * j * (2n + 2, n + 1 - j ) / (2n, n)

    //FUDO| we start with counting with j=0, but the coefficients are B_(j+1)

    int ordpo = len + 1;

    //FUDO| double-check all the coefficients
    for ( int i=0; i<nord; i++ ) {
        int k = i+1;
        coeffs[i] = get_Bj (ordpo, k);
    }
}

double FixASPC::get_Bj(int n, int k)
{
    return pow(-1, k+1) * k * combi (2 * n + 2, n + 1 - k ) / combi (2 * n, n);
}

// FUDO| taken from: http://stackoverflow.com/questions/24294192/computing-the-binomial-coefficient-in-c
long long FixASPC::combi(int n,int k)
{
    long long ans=1;
    k=k>n-k?n-k:k;
    int j=1;
    for(;j<=k;j++,n--)
    {
        if(n%j==0)
        {
            ans*=n/j;
        }else
        if(ans%j==0)
        {
            ans=ans/j*n;
        }else
        {
            ans=(ans*n)/j;
        }
    }
    return ans;
}

/* ----------------------------------------------------------------------
   allocate/initialize memory for a new vector with N elements per atom
------------------------------------------------------------------------- */

void FixASPC::add_vector(int n)
{
  memory->grow(peratom,nvector+1,"aspc:peratom");
  peratom[nvector] = n;

  vectors = (double **)
    memory->srealloc(vectors,(nvector+1)*sizeof(double *),"aspc:vectors");

  coeffs = (double *)
    memory->srealloc(coeffs,(nvector+1)*sizeof(double),"aspc:coeffs");

  memory->create(vectors[nvector],atom->nmax*n,"aspc:vector");

  int ntotal = n*atom->nlocal;
  for (int i = 0; i < ntotal; i++) vectors[nvector][i] = 0.0;
  nvector++;
}

/* ----------------------------------------------------------------------
   return a pointer to the Mth vector
------------------------------------------------------------------------- */

double *FixASPC::request_vector(int m)
{
  return vectors[m];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixASPC::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvector; m++)
    bytes += atom->nmax*peratom[m]*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixASPC::grow_arrays(int nmax)
{
  for (int m = 0; m < nvector; m++)
    memory->grow(vectors[m],peratom[m]*nmax,"aspc:vector");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixASPC::copy_arrays(int i, int j, int delflag)
{
  int m,iper,nper,ni,nj;

  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*i;
    nj = nper*j;
    for (iper = 0; iper < nper; iper++) vectors[m][nj++] = vectors[m][ni++];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixASPC::pack_exchange(int i, double *buf)
{
  int m,iper,nper,ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*i;
    for (iper = 0; iper < nper; iper++) buf[n++] = vectors[m][ni++];
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixASPC::unpack_exchange(int nlocal, double *buf)
{
  int m,iper,nper,ni;

  int n = 0;
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    ni = nper*nlocal;
    for (iper = 0; iper < nper; iper++) vectors[m][ni++] = buf[n++];
  }
  return n;
}

//FUDO| AFAICT we also need pack_forward and unpack_forward?!?
