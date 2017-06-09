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
#include "fix_aspc_drude.h"
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

enum{MAXITER,MAXEVAL,FTOL};

/* ---------------------------------------------------------------------- */

FixASPCDrude::FixASPCDrude(LAMMPS *lmp, int narg, char **arg) : FixASPC(lmp,narg,arg)
{

  if (narg < 7) error->all(FLERR,"Illegal fix ASPCDrude command");

  kd = force->numeric(FLERR,arg[5]);

  int n = strlen(arg[6]) + 1;
  id_ef = new char[n];
  strcpy(id_ef,arg[6]);

  what = COORDS;
  dim = 3;

  neval = 1;
  scf = 0;
  ftol = 1.e-02;
  printconv = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scf") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix aspc/drude command");
      scf = 1;
      ftol = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"neval") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal  aspc/drude command");
      neval = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix aspc/drude command");
  }

  int natom = atom->natoms;
  memory->create(dx, (natom+1)*dim, "aspc_drude:dx");
  memory->create(hrm, (natom+1)*dim, "aspc_drude:hrm");

}

/* ---------------------------------------------------------------------- */

FixASPCDrude::~FixASPCDrude()
{
  delete [] id_ef;

  memory->destroy(dx);
  memory->destroy(hrm);
}

void FixASPCDrude::init()
{
  FixASPC::init();

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"drude") == 0) break;
  if (ifix == modify->nfix) error->all(FLERR, "ASPCDrude for DRUDE requires fix drude");
  fix_drude = (FixDrude *) modify->fix[ifix];

  //FUX| handle if compute id not found
  int ief = modify->find_compute(id_ef);
  if (ief < 0)
    error->all(FLERR,"Compute ID for fix efield/atom does not exist");

  c_ef = modify->compute[ief];

  //FU| perform additional check whether only drude types have been included

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  int *drudetype = fix_drude->drudetype;

  for ( int i=0; i<nlocal; i++ )
    if ( mask[i] & groupbit && (drudetype[type[i]] != DRUDE_TYPE) )
      error->all(FLERR, "Your fix ASPC group contains non-Drude particles");

}

void FixASPCDrude::setup_pre_force(int vflag)
{
    //FU| this needs to be included, and somewhere else than in class construction, or init()
    drudeid = fix_drude->drudeid;

    FixASPC::setup_pre_force(vflag);
}

void FixASPCDrude::correct()
{
    int i, k, n, baseind;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;

    int coreind;
    double df;

    double onemdamp = 1. - damp;

    c_ef->compute_peratom();
    
    //FUX | AFAICT the array shouldn't grow in one correct() step
    double **f = c_ef->array_atom;

    baseind = 0;

    r2r_indices_forward();

    //FU| save the predicted positions
    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {

        coreind = (int) drudeid[i];
        coreind *= dim;

        for ( k=0; k<dim; k++)
            dx[baseind+k] = qty[baseind+k] - qty[coreind+k];

      }

      baseind += dim;
    }

    baseind = 0;
    int conv = 0;
    int nmaxener = neval - 1;

    //FU| apply corrector neval times
    for ( n=0; n<neval; n++ ) {
        baseind = 0;

        for ( i=0; i<nlocal; i++ ) {
          if ( mask[i] & groupbit ) {

            coreind = (int) drudeid[i];
            coreind *= dim;

            for ( k=0; k<dim; k++) {

              //FUX| calculate new dipole moments, assumes that fieldforce option is set for compute efield/atom
              df = f[i][k] / kd;
              qty[baseind+k] = qty[coreind+k] + df;
            }

          }

          baseind += dim;
        }

        //FU| additional communication if multiple corrector applications, otherwise one in pre_force, after correct()
        if ( (scf) || ( neval > 1) )
          comm->forward_comm();

        if ( scf ) {
          r2r_indices_reverse();
          c_ef->compute_peratom();
          r2r_indices_forward();

          //FU| we don't need the energy, convergence only checked on forces
          calc_spring_forces_energy();

          conv = check_convergence(f);

          //FU| possibly add other convergence criteria
          if ( conv ) {
            if ( ( printconv ) && ( comm->me == 0) ) {
              if (conv == FTOL) printf("Forces ");

              printf("converged in %i iterations\n", n+1);
            }
            break;
          }

        }
        else if ( n < nmaxener ) {
          r2r_indices_reverse();
          c_ef->compute_peratom();
          r2r_indices_forward();
        }

    }

    if ((scf) && !(conv))
      error->all(FLERR,"Convergence could not be achieved in maximum number of iterations");

    baseind = 0;

    //FU| perform the actual ASPC step, between predictor and corrector
    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {

        coreind = (int) drudeid[i];
        coreind *= dim;

        for ( k=0; k<dim; k++) {
            df = qty[baseind+k] - qty[coreind+k];
            qty[baseind+k] = qty[coreind+k] + damp * df + onemdamp * dx[baseind+k];
        }
      }
      baseind += dim;
    }

    force_clear();

    r2r_indices_reverse();
}

void FixASPCDrude::r2r_indices_forward()
{
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i=0; i<nlocal; i++)
      if (mask[i] & groupbit )
        drudeid[i] = (tagint) domain->closest_image(i, atom->map(drudeid[i]));
}

void FixASPCDrude::r2r_indices_reverse()
{
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i=0; i<nlocal; i++)
      if (mask[i] & groupbit)
        drudeid[i] = atom->tag[(int) drudeid[i]];
}

//FU| additional prediction routine, because we need the position difference and not the absolute position
//FUDO| this could be shortened by calling  ASPC::predict() and then adding the corresponding values to the core positions
void FixASPCDrude::predict()
{

    int i, k, baseind, coreind;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;

    comm->forward_comm();
    r2r_indices_forward();

    //FU| use the general ASPC prediction
    FixASPC::predict();

    //FU| and add the core position
    //FUDO| additional communication needed?

    baseind = 0;

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {
        coreind = (int) drudeid[i];
        coreind *= dim;

        for ( k=0; k<dim; k++)
           qty[baseind+k] += qty[coreind+k];

      }

      baseind += dim;
    }

    r2r_indices_reverse();

}

void FixASPCDrude::cpy2hist()
{
    int i, k, coreind, baseind;
    int nlocal = atom->nlocal;
    double *data;
    int *mask = atom->mask;

    data = request_vector(0);

    r2r_indices_forward();

    baseind = 0;

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {
        coreind = (int) drudeid[i];
        coreind *= dim;

        for ( k=0; k<dim; k++ )
          data[baseind+k] = qty[baseind+k] - qty[coreind+k];
      }
      baseind += dim;
    }

    r2r_indices_reverse();

}

int FixASPCDrude::check_convergence(double **f)
{
  int nlocal = atom->nlocal;
  int i;
  int *mask = atom->mask;
  int noconv = 0;
  int allnoconv = 0;
  int baseind;

  double fsqr;
  double ftolsqr = ftol*ftol;
  double frc[dim];

  baseind = 0;
  for ( i=0; i<nlocal; i++ ) {
    if ( mask[i] & groupbit ) {
  
      // forces are assumed, i.e., fieldforce option to compute efield/atom is on
      frc[0] = f[i][0] + hrm[baseind];
      frc[1] = f[i][1] + hrm[baseind+1];
      frc[2] = f[i][2] + hrm[baseind+2];

      fsqr = frc[0]*frc[0] + frc[1]*frc[1] + frc[2]*frc[2];

      if (fsqr > ftolsqr)
        noconv += 1;
    }
    baseind += dim;
  }

  MPI_Allreduce(&noconv,&allnoconv,1,MPI_INT,MPI_SUM,universe->uworld);

  if ( !(allnoconv) )
    return FTOL;

  return 0;

}

double FixASPCDrude::calc_spring_forces_energy()
{
  int i, coreind, baseind;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delx, dely, delz, fbond, ebond;
  double rsq,r,dr,rk;

  baseind = 0;
  ebond = 0.;

  for ( i=0; i<nlocal; i++ ) {
    if ( mask[i] & groupbit ) {
  
      coreind = (int) drudeid[i];
      coreind *= dim;

      delx = qty[baseind+0] - qty[coreind+0];
      dely = qty[baseind+1] - qty[coreind+1];
      delz = qty[baseind+2] - qty[coreind+2];

      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);
      dr = r;
      rk = kd * dr;

      // force & energy

      ebond += rk*dr;
      fbond = -rk/r;

      hrm[baseind] = delx * fbond;
      hrm[baseind+1] = dely * fbond;
      hrm[baseind+2] = delz * fbond;
    }
    baseind += dim;
  }

  return ebond;

}

int FixASPCDrude::modify_param(int narg, char **arg)
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
