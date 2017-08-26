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
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_aspc.h"
#include "fix_aspc_iccs.h"
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
#include "universe.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPS_ENERGY 1.0e-8
#define TWOPI 6.283185307179586
#define EPS 1.E-02

enum{MAXITER,MAXEVAL,FTOL};

/* ---------------------------------------------------------------------- */

FixASPCICCS::FixASPCICCS(LAMMPS *lmp, int narg, char **arg) : FixASPC(lmp,narg,arg)
{

  if (narg < 12) error->all(FLERR,"Illegal fix ASPCICCS command");

  //FUX| need input:
  //      - compute/efield
  //      - property names: 1 dielectric 3 components surface normal

  bulk_perm = force->numeric(FLERR,arg[5]);

  int n = strlen(arg[6]) + 1;
  id_ef = new char[n];
  strcpy(id_ef,arg[6]);

  n = strlen(arg[7]) + 1;
  id_diel = new char[n];
  strcpy(id_diel, arg[7]);

  n = strlen(arg[8]) + 1;
  id_area = new char[n];
  strcpy(id_area, arg[8]);

  n = strlen(arg[9]) + 1;
  id_srfx = new char[n];
  strcpy(id_srfx, arg[9]);

  n = strlen(arg[10]) + 1;
  id_srfy = new char[n];
  strcpy(id_srfy, arg[10]);

  n = strlen(arg[11]) + 1;
  id_srfz = new char[n];
  strcpy(id_srfz, arg[11]);

  comm_forward = 1;

  what = CHARGES;
  dim = 1;

  neval = 1;
  scf = 0;
  nfail = 0;
  faild = 0;
  printconv = 0;
  qinit = 0;
  conv = EPS;
  recalcf = 0;
  reinitialize = 1;

  int iarg = 12;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scf") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      scf = 1;
      conv = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"neval") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      neval = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

}

/* ---------------------------------------------------------------------- */

FixASPCICCS::~FixASPCICCS()
{
  delete [] id_ef;
  delete [] id_diel;
  delete [] id_area;
  delete [] id_srfx;
  delete [] id_srfy;
  delete [] id_srfz;

  memory->destroy(qprv);
  memory->destroy(qprd);
  memory->destroy(qnxt);
}

void FixASPCICCS::init()
{
  FixASPC::init();

  //FUX| handle if compute id not found
  int ief = modify->find_compute(id_ef);
  if (ief < 0)
    error->all(FLERR,"Compute ID for fix efield/atom does not exist");

  c_ef = modify->compute[ief];

  int natoms = atom->natoms;

  qprv = memory->create(qprv,natoms+1,"iccs:qprv");
  qprd = memory->create(qprd,natoms+1,"iccs:qprd");
  qnxt = memory->create(qnxt,natoms+1,"iccs:qnxt");

  add_vector(dim);
  contrast = request_vector(nvector-1);

  reset_vectors();
  calculate_contrast();
}

void FixASPCICCS::setup_pre_force(int vflag)
{
    //FU| this needs to be included, and somewhere else than in class construction, or init()

  reset_vectors();
  initialize_charges();
  FixASPC::setup_pre_force(vflag);
}

void FixASPCICCS::reset_vectors()
{
  // necessary?
  // FixASPC::reset_vectors();

  int index, flag;
  
  index = atom->find_custom(id_diel, flag);
  p_diel = atom->dvector[index];
  
  index = atom->find_custom(id_area, flag);
  p_area = atom->dvector[index];
  
  index = atom->find_custom(id_srfx, flag);
  p_srfx = atom->dvector[index];
  
  index = atom->find_custom(id_srfy, flag);
  p_srfy = atom->dvector[index];

  index = atom->find_custom(id_srfz, flag);
  p_srfz = atom->dvector[index];
  contrast = request_vector(nvector-1);
}

void FixASPCICCS::correct()
{
  int n;

  // c_ef->compute_peratom();

  //FUX | AFAICT the array shouldn't grow in one correct() step
  // double **f = c_ef->array_atom;

  //FU| apply corrector neval times

  if ( (recalcf) || !(c_ef->invoked_peratom == update->ntimestep) )
      c_ef->compute_peratom();

  reset_vectors();

  int converged = 0;

  int halfeval = neval / 2;

  double dampiter = damp / 11.;
  double dampbak = damp;

  int tentheval = neval / 10;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *q = atom->q;

  for( int i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      qprd[i] = q[i];

RESTARTDAMP:
  for ( n=0; n<neval; n++ ) {

    // includes charge update and communication
    // printf("HELLO THERE\n");
    iterate();

    if ( scf ) {
      converged = check_convergence();

      if( converged == 2 )
        if( comm->me == 0 )
          printf("Some charges are suspiciously close to zero\n");

      //FU| possibly add other convergence criteria
      if ( converged ) {
        if ( ( printconv ) && ( comm->me == 0) ) {
          printf("converged in %i iterations\n", n+1);
        }
        break;
      }

      c_ef->compute_peratom();

    }

    // if ( !(converged) && (n == halfeval) && (reinitialize) ) {
    //   if ( comm->me == 0 )
    //     printf("Re-initialization charges on ICC* particles\n");

    //   initialize_charges();
    // }

    // printf("%i %i\n", n % 10, tentheval);
    if( (scf) && !(converged) && (n == neval-1) && (reinitialize) )
    {
      //FUX| TODO need to check that damping isn't too low, otherwise there will be no measurable change in charges
      damp -= dampiter;
      printf("Resetting damping factor to %g\n", damp);
      // for( int i=0; i<nlocal; ++i )
      //   if( mask[i] & groupbit )
      //     q[i] = qprd[i];

      // comm->forward_comm_fix(this);
      // force->kspace->qsum_qsq();

      // this would break time-reversibility?
      initialize_charges();
      goto RESTARTDAMP;
    }

  }

  damp = dampbak;

  if ((scf) && !(converged) && (faild >= nfail))
    error->all(FLERR,"Convergence could not be achieved in maximum number of iterations");
  else {
    faild += 1;
  }

  // checkme();

}

void FixASPCICCS::checkme()
{
  int *mask = atom->mask;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  double qtot = 0.;
  double allqtot = 0.;

  for( int i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      qtot += q[i];

  MPI_Allreduce(&qtot,&allqtot,1,MPI_DOUBLE,MPI_SUM,world);
  if( comm->me == 0 ) printf("My total ICC* charge is: %g\n", qtot);
}

void FixASPCICCS::iterate()
{
  backup_charges();

  //FUX| this could become more complicated if fixes like ASPC/DRUDE SCF calculations get involved
  calculate_charges_iccs();
  update_charges();
}

void FixASPCICCS::backup_charges()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double *q = atom->q;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      qprv[i] = q[i];
}

//FUX| calculate_charges()
//FUX| update_charges()

void FixASPCICCS::calculate_charges_iccs()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  double **ef = c_ef->array_atom;

  // printf("NEW ITERATION\n");
  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      qnxt[i] = contrast[i] * ( ef[i][0]*p_srfx[i] + ef[i][1]*p_srfy[i] + ef[i][2]*p_srfz[i] );

}

void FixASPCICCS::initialize_charges()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      q[i] = 0.01 * ( (float) rand() / RAND_MAX - 0.5);

  comm->forward_comm_fix(this);
  if( kspace_compute_flag )
    force->kspace->qsum_qsq();
      
  qinit = 1;
}

void FixASPCICCS::update_charges()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;
  double onemdamp = 1. - damp;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit )
      q[i] = onemdamp * qprv[i] + damp * qnxt[i];

  comm->forward_comm_fix(this);
  if( kspace_compute_flag )
    force->kspace->qsum_qsq();
}

void FixASPCICCS::predict()
{
  FixASPC::predict();
  comm->forward_comm_fix(this);
  if( kspace_compute_flag )
    force->kspace->qsum_qsq();
}

int FixASPCICCS::check_convergence()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *q = atom->q;

  int isnotconv = 0;
  int allisnotconv = 0;

  int isdangerous = 0;
  int allisdangerous = 0;

  for( i=0; i<nlocal; ++i )
    if( mask[i] & groupbit ) {
      if( fabs( ( q[i] - qprv[i] ) / qprv[i] ) > conv )
        isnotconv += 1;

      if( fabs( q[i] ) < EPS_ENERGY )
        isdangerous += 1;
    }

  //FUX| substitute MPI_SUM by MPI_MAX? and check for largest charge
  //FUX| avoids havhing to calculate the sum...
  MPI_Allreduce(&isnotconv,&allisnotconv,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&isdangerous,&allisdangerous,1,MPI_INT,MPI_SUM,world);

  // if( comm->me == 0 ) printf("ICCS noconv is: %i\n", allisnotconv);

  // if( comm->me == 0 ) printf("ALLISNOTCONV: %i\n", allisnotconv);
  if ( allisnotconv )
    return 0;
  else {
    if ( allisdangerous )
      return 2;
    else
      return 1;
  }
}

void FixASPCICCS::calculate_contrast()
{
  int i;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  // double fpieps = 0.0030119505336064496;
  double fpieps = 0.0030114702717118813;

  for ( i=0; i<nlocal; i++ )
    if( mask[i] & groupbit ) {

      contrast[i] = bulk_perm / TWOPI * p_area[i] * fpieps;

      if ( p_diel[i] < 1 )
        contrast[i] *= -1.;
      else
        contrast[i] *= (bulk_perm - p_diel[i]) / (bulk_perm + p_diel[i]);
    }

}

int FixASPCICCS::modify_param(int narg, char **arg)
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
    } else if (strcmp(arg[iarg],"neval") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      neval = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"printconv") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1], "yes") == 0) printconv = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) printconv = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"scf") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1], "yes") == 0) scf = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) scf = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"failures") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      nfail = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"conv") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      conv = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"recalcf") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1], "yes") == 0) recalcf = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) recalcf = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"reinitialize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      if (strcmp(arg[iarg+1], "yes") == 0) reinitialize = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) reinitialize = 0;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

  return 2;
}

/* ---------------------------------------------------------------------- */

int FixASPCICCS::pack_forward_comm(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int m;

  for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];

  return n;
}

/* ---------------------------------------------------------------------- */

void FixASPCICCS::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  for(m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}
