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
#include "fix_aspc_dipole.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_poldip.h"
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

enum{MAXITER,MAXEVAL,FTOL};

/* ---------------------------------------------------------------------- */

FixASPCDipole::FixASPCDipole(LAMMPS *lmp, int narg, char **arg) : FixASPC(lmp,narg,arg)
{

  if (narg < 5) error->all(FLERR,"Illegal fix ASPCDipole command");

  what = DIPOLES;
  dim = 4;
  ndx_dim = 3;

  // FUDO| need additional check whether lj/cut/coul/<cut/long>/sa is active, if not abort...
  // make sure dim is always the right number here!

  neval = 1;
  scf = 0;
  ftol = 1.e-02;
  printconv = 0;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scf") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      scf = 1;
      ftol = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"neval") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      neval = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

  int natom = atom->natoms;
  memory->create(dx, (natom+1)*dim, "aspc_dipole:dx");
  memory->create(oldf, (natom+1)*dim, "aspc_dipole:oldf");

  fpieps = // 0.0030116665;
           // 0.0030119505336064496;
           0.0030114702717118813;
           /* convert 4 * pi * eps_0 (farad/m, coulomb/volt/m, coulomb^2/joule/m) to kcal/mol/e/angstrom (this is the cgs conversion factor) */

}

/* ---------------------------------------------------------------------- */

FixASPCDipole::~FixASPCDipole()
{
  memory->destroy(dx);
  memory->destroy(oldf);
}

void FixASPCDipole::init()
{
  FixASPC::init();

  int ntotal = dim*atom->natoms;
  for (int i = 0; i < ntotal; i++) oldf[i] = 0.0;
}

void FixASPCDipole::correct()
{
    int i, k, n, baseind;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;

    double df;

    double eprevious, ecurrent;

    double **f = atom->f;
    double onemdamp = 1. - damp;
    double *alf = atom->alf;

    eprevious = energy_force_es(1);

    baseind = 0;

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {

        //FU| save the predicted dipoles
        for ( k=0; k<ndx_dim; k++)
            dx[baseind+k] = qty[baseind+k];

      }

      baseind += dim;
    }

    baseind = 0;
    int conv = 0;
    int nmaxener = neval - 1;

    for ( n=0; n<neval; n++ ) {
        baseind = 0;

        double nrm;
        for ( i=0; i<nlocal; i++ ) {
          if ( mask[i] & groupbit ) {

          //FU| this might lead to issues when initially a point dipole is defined for the atom, but no polarizability
          //FU| that is a case we cannot treat right now
          if ( alf[i] == 0.0 ) {
              baseind += dim;
            continue;
          }

            nrm = 0.;
            for ( k=0; k<ndx_dim; k++) {

              df = fpieps * alf[i] * f[i][k];  // q[i];
              qty[baseind+k] = df;
              nrm += df*df;
            }
            nrm = sqrt(nrm);
            qty[baseind+k] = nrm;

          }

          baseind += dim;
        }

        if ( (scf) || ( neval > 1) )
          comm->forward_comm();

        if ( scf ) {
          ecurrent = energy_force_es(1);

          conv = check_convergence(ecurrent, eprevious, f, oldf);
          eprevious = ecurrent;

          if ( conv ) {
            if ( ( printconv ) && ( comm->me == 0) ) {
              if (conv == FTOL) printf("Forces ");

              printf("converged in %i iterations\n", n+1);
            }
            break;
          }

        }
        else if ( n < nmaxener ) {
          ecurrent = energy_force_es(1);
        }

    }

    if ((scf) && !(conv))
      error->all(FLERR,"Convergence could not be achieved in maximum number of iterations");

    baseind = 0;
    double nrm;

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {
        nrm = 0.;
        for ( k=0; k<ndx_dim; k++) {
            qty[baseind+k] = damp * qty[baseind+k] + onemdamp * dx[baseind+k];
            nrm += qty[baseind+k]*qty[baseind+k];
        }

        //FU| update the total dipole moment
        //FU| it should be k here, because it already has been increased, should it not?
        qty[baseind+k] = sqrt(nrm);

      }
      baseind += dim;
    }

    force_clear();

}

//FU| this could be adapted to new variables, i.e., check change in dipole moment
int FixASPCDipole::check_convergence(double ecurrent, double eprevious, double **f, double *oldf)
{
  int nlocal = atom->nlocal;
  int i;
  int *mask = atom->mask;
  int noconv = 0;
  int allnoconv = 0;
  int baseind;

  double fsqr;
  double ftolsqr = ftol*ftol;
  double dltf[ndx_dim];

  baseind = 0;
  for ( i=0; i<nlocal; i++ ) {
    if ( mask[i] & groupbit ) {
  
      dltf[0] = f[i][0] - oldf[baseind+0];
      dltf[1] = f[i][1] - oldf[baseind+1];
      dltf[2] = f[i][2] - oldf[baseind+2];

      fsqr = dltf[0]*dltf[0] + dltf[1]*dltf[1] + dltf[2]*dltf[2];

      if (fsqr > ftolsqr)
        noconv += 1;

      //FU| update old field
      oldf[baseind+0] = f[i][0];
      oldf[baseind+1] = f[i][1];
      oldf[baseind+2] = f[i][2];

    }
    baseind += dim;
  }

  MPI_Allreduce(&noconv,&allnoconv,1,MPI_INT,MPI_SUM,universe->uworld);

  if ( !(allnoconv) )
    return FTOL;

  return 0;

}

double FixASPCDipole::energy_force_es(int resetflag)
{
  force_clear();

  timer->stamp();

  //FU| make sure to only get electric field

  force->pair->compute_efld = 1;
  if (pair_compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }

  if (kspace_compute_flag) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }

  if (modify->n_pre_reverse) {
    modify->pre_reverse(eflag,vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
  }

  force->pair->compute_efld = 0;

  double energy = compute_pe_scalar();

  return energy;
}

double FixASPCDipole::compute_pe_scalar()
{
  double scalar = 0.;

  double one = 0.0;
  if (force->pair)
    one += force->pair->eng_vdwl + force->pair->eng_coul;

  if (atom->molecular) {
    if (force->bond) one += force->bond->energy;
    if (force->angle) one += force->angle->energy;
    if (force->dihedral) one += force->dihedral->energy;
    if (force->improper) one += force->improper->energy;
  }

  MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->kspace) scalar += force->kspace->energy;

  if (force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    scalar += force->pair->etail / volume;
  }

  return scalar;
}

int FixASPCDipole::modify_param(int narg, char **arg)
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
    } else if (strcmp(arg[iarg],"ftol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix_modify command");
      ftol = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_modify command");
  }

  return 2;
}

void FixASPCDipole::predict()
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
    }
    else
        my_nord = nord;

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {
        for ( k=0; k<ndx_dim; k++) {
           qty[baseind+k] = coeffs[0] * data[baseind+k];
        }
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

    //FU| calculating the total dipole moment seperatly, not predicted!
    double nrm;
    baseind = 0;
    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {
        nrm = 0.;
        for ( k=0; k<ndx_dim; k++) {
           nrm += qty[baseind+k]*qty[baseind+k];
        }
        nrm = sqrt(nrm);
        qty[baseind+k] = nrm;
      }

      baseind += dim;
    }
}

//FUDO| need function that adds polarization energy to the total energy, if requested (i.e., 0.5 mu^2/alpha)
