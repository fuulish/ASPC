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

  //FU| defaults
  what = NONE;
  dim = 0;

  eflag = 1;
  vflag = 0;

  first = update->ntimestep;
  nonlyhist = 0;
  filledhist = 0;

  //FU| memory-related
  nvector = 0;
  peratom = NULL;
  vectors = NULL;
  coeffs = NULL;

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
   // FU| this should be included in some way, not sure how
   
   if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
   else pair_compute_flag = 0;
   if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
   else kspace_compute_flag = 0;

   bigint ndofme = 3 * static_cast<bigint>(atom->nlocal);
   MPI_Allreduce(&ndofme,&ndoftotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

   if ( what == NONE )
      error->all(FLERR, "Don't know what to predict/correct");

   for (int k=0; k<nord; k++)
     add_vector(dim);

   tstart = update->ntimestep;

   generate_coefficients(length);

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
/* FU| these could/should be included in later versions
            case DIPOLE:
                qty = atom->mu;
                break;
            case QEQ:
                qty = atom->q;
                break;
            case QEQREAX:
                qty = atom->q;
                break;
*/
            default:
                error->all(FLERR, "Do not know which quantity to predict/correct");
                break;
        }
    }

}

void FixASPC::update_history()
{
    //FU| move current state to history and all other previous histories

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
    //FU| copy current state to history

    int i;
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
    if (strcmp(arg[iarg],"kspace_compute_flag") == 0) {
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

void FixASPC::final_integrate()
{
    reset_vectors();
    update_history();
}

void FixASPC::initial_integrate(int vflag)
{

    reset_vectors();

    if ( !(filledhist) ) {
      update_history();

      filledhist = 1;
    }

    if ( ( update->ntimestep - tstart ) < (nonlyhist))
        return;

    predict();
}

void FixASPC::setup_pre_force(int vflag)
{
    //FU| there's nothing in the history yet, and nothing to predict from, hence only to the correction
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

    //FU| we'll not do re-neighboring, should we? 
    comm->forward_comm();
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

    //FU| if we don't have a fully accumulated history, yet, start with a smaller one
    if ( tlength <= nord ) {
        my_nord = tlength;
        generate_coefficients(my_nord - 2);
    }
    else
        my_nord = nord;

    //FU| use the first iteration as assignment, would also work if qty gets reset to zero each time
    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit )
        for ( k=0; k<dim; k++) {
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

//FU| specific routine needs to be added for each fix
void FixASPC::correct()
{
}

int FixASPC::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

void FixASPC::generate_coefficients(int len)
{
    //FU| we start with counting with j=0, but the coefficients are B_(j+1)

    int ordpo = len + 1;

    for ( int i=0; i<nord; i++ ) {
        int k = i+1;
        coeffs[i] = get_Bj (ordpo, k);
    }
}

double FixASPC::get_Bj(int n, int k)
{
    // B_j = (-1)^(j+1) * j * (2n + 2, n + 1 - j ) / (2n, n)

    return pow(-1, k+1) * k * combi (2 * n + 2, n + 1 - k ) / combi (2 * n, n);
}

// FU| taken from: http://stackoverflow.com/questions/24294192/computing-the-binomial-coefficient-in-c
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
