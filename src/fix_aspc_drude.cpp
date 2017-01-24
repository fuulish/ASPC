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

  if (narg < 6) error->all(FLERR,"Illegal fix ASPCDrude command");

  printf("%i number of arguments\n", narg);

  kd = force->numeric(FLERR,arg[5]);
  printf("HARMONIC FORCE CONSTANT IS: %14.8f\n", kd);
  // FUDO| check which way we wanna handle things here, either LAMMPS-equivalent of force constant, or 1./2. * k
  // FUDO| possible conversion factors?!
  // kd *= 2.;

  what = COORDS;
  dim = 3;

  // FUDO| need additional check whether lj/cut/coul/<cut/long>/sa is active, if not abort...
  // make sure dim is always the right number here!

  zerovels = 1;
  neval = 1;
  scf = 0;
  ftol = 1.e-02;
  printconv = 0;

  int natom = atom->natoms;
  memory->create(dx, (natom+1)*dim, "aspc_drude:dx");
  memory->create(hrm, (natom+1)*dim, "aspc_drude:hrm");

  // comm_forward = 9;
}

/* ---------------------------------------------------------------------- */

FixASPCDrude::~FixASPCDrude()
{
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

  //FUDO| perform additional check whether only drude types have been included

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  int *drudetype = fix_drude->drudetype;

  // FUDO| alternatively loop only over drude atoms in predict() and correct(), opinions?
  for ( int i=0; i<nlocal; i++ )
    if ( mask[i] & groupbit && (drudetype[type[i]] != DRUDE_TYPE) )
      error->all(FLERR, "Your fix ASPC group contains non-Drude particles");

}

void FixASPCDrude::setup_pre_force(int vflag)
{
    //FUDO| only works if drudeid is set here, or in correct/predict-kinda location
    reset_vectors();
    drudeid = fix_drude->drudeid;

    // pre_force(vflag);
    // comm->forward_comm();

    update_history();

}

void FixASPCDrude::correct()
{
    int i, k, n, baseind;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;

    int coreind;
    double df;

    double eprevious, ecurrent;

    double **f = atom->f;
    double *q = atom->q;
    double onemdamp = 1. - damp;

    // double **x = atom->x;
    // for ( int l=0; l<nlocal; l++ )
    //   if ( mask[l] & groupbit )
    //     printf("%14.8f %14.8f %14.8f\n", x[l][0], x[l][1], x[l][2]);
    //     // printf("%14.8f %14.8f %14.8f\n", f[l][0], f[l][1], f[l][2]);

    eprevious = energy_force_es(1);
    // eprevious += calc_spring_forces_energy();

    // FUDO| maybe remove re-neighboring from energy_force_es
    // printf("BEFORE ENERGY\n");
    // printf("AFTER ENERGY\n");

    baseind = 0;

    r2r_indices_forward();

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {

        coreind = (int) drudeid[i];
        coreind *= dim;

        // FUDO| save the predicted positions
        for ( k=0; k<dim; k++)
            dx[baseind+k] = qty[baseind+k] - qty[coreind+k];

      }

      baseind += dim;
    }

    baseind = 0;
    int conv = 0;
    int nmaxener = neval - 1;

    for ( n=0; n<neval; n++ ) {
        baseind = 0;

        for ( i=0; i<nlocal; i++ ) {
          if ( mask[i] & groupbit ) {

            coreind = (int) drudeid[i];
            coreind *= dim;

            for ( k=0; k<dim; k++) {

              df = f[i][k] / kd;
              qty[baseind+k] = qty[coreind+k] + df;
            }

          }

          baseind += dim;
        }

        if ( (scf) || ( neval > 1) )
          comm->forward_comm();

        if ( scf ) {
          //FUDO| also calculate energy contribution of harmonic springs?
          r2r_indices_reverse();
          ecurrent = energy_force_es(1);
          r2r_indices_forward();

          calc_spring_forces_energy();
          // ecurrent += calc_spring_forces_energy();

          conv = check_convergence(ecurrent, eprevious, f);
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
          r2r_indices_reverse();
          ecurrent = energy_force_es(1);
          r2r_indices_forward();
        }

    }

    if ((scf) && !(conv))
      error->all(FLERR,"Convergence could not be achieved in maximum number of iterations");

    baseind = 0;

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit ) {

        coreind = (int) drudeid[i];
        coreind *= dim;

        //FUDO| here we again need to make sure that we're using the difference between core and dp

        for ( k=0; k<dim; k++) {
            df = qty[baseind+k] - qty[coreind+k];
            qty[baseind+k] = qty[coreind+k] + damp * df + onemdamp * dx[baseind+k];
        }
        // double **x = atom->x;
        // printf("%5i %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n", i, f[i][0], f[i][1], f[i][2], x[i][0], x[i][1], x[i][2]);
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

void FixASPCDrude::predict()
{

    int i, k, n, baseind, coreind;
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
        // printf("TEMPORARY: %i at TIMESTEP: %5i\n", tlength, update->ntimestep);
        // for ( n=0; n < my_nord; n++ )
        //     printf("COEFF: %14.8f\n", coeffs[n]);
    }
    else
        my_nord = nord;

    // comm->forward_comm_fix(this);

    comm->forward_comm();
    r2r_indices_forward();

    for ( i=0; i<nlocal; i++ ) {
      if ( mask[i] & groupbit )
        for ( k=0; k<dim; k++)
           qty[baseind+k] = coeffs[0] * data[baseind+k];

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

    // comm->forward_comm_fix(this);
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

int FixASPCDrude::check_convergence(double ecurrent, double eprevious, double **f)
{
  int nlocal = atom->nlocal;
  int i, k;
  int *mask = atom->mask;
  int noconv = 0;
  int allnoconv = 0;
  int coreind, baseind;

  double fsqr;
  double ftolsqr = ftol*ftol;
  double frc[dim];

  baseind = 0;
  for ( i=0; i<nlocal; i++ ) {
    if ( mask[i] & groupbit ) {
  
      frc[0] = f[i][0] + hrm[baseind];
      frc[1] = f[i][1] + hrm[baseind+1];
      frc[2] = f[i][2] + hrm[baseind+2];

      fsqr = frc[0]*frc[0] + frc[1]*frc[1] + frc[2]*frc[2];
      // printf("%14.8f %14.8f\n", f[i][0], delx*fbond);

      if (fsqr > ftolsqr)
        noconv += 1;
    }
    baseind += dim;
  }

  // printf("%14.8f %14.8f %i\n", ecurrent, eprevious, noconv);

  MPI_Allreduce(&noconv,&allnoconv,1,MPI_INT,MPI_SUM,universe->uworld);

  if ( !(allnoconv) )
    return FTOL;

  return 0;

}

double FixASPCDrude::energy_force_es(int resetflag)
{
  force_clear();

  timer->stamp();

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

  double energy = compute_pe_scalar();

  return energy;
}

double FixASPCDrude::compute_pe_scalar()
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

/* ---------------------------------------------------------------------- */
int FixASPCDrude::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  double ** x = atom->x, ** v = atom->v, ** f = atom->f;
  int * type = atom->type, * drudetype = fix_drude->drudetype;
  double dx,dy,dz;
  int dim = domain->dimension;
  int m = 0;
  for (int i=0; i<n; i++) {
    int j = list[i];
    if (pbc_flag == 0 ||
        (fix_drude->is_reduced && drudetype[type[j]] == DRUDE_TYPE)) {
        for (int k=0; k<dim; k++) buf[m++] = x[j][k];
    }
    else {
        if (domain->triclinic != 0) {
            dx = pbc[0]*domain->xprd + pbc[5]*domain->xy;
            dy = pbc[1]*domain->yprd;
            if (dim == 3) {
                dx += + pbc[4]*domain->xz;
                dy += pbc[3]*domain->yz;
                dz = pbc[2]*domain->zprd;
            }
        }
        else {
            dx = pbc[0]*domain->xprd;
            dy = pbc[1]*domain->yprd;
            if (dim == 3)
                dz = pbc[2]*domain->zprd;
        }
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        if (dim == 3)
            buf[m++] = x[j][2] + dz;
    }
    for (int k=0; k<dim; k++) buf[m++] = v[j][k];
    for (int k=0; k<dim; k++) buf[m++] = f[j][k];
  }
  return m;
}

/* ---------------------------------------------------------------------- */
void FixASPCDrude::unpack_forward_comm(int n, int first, double *buf)
{
  double ** x = atom->x, ** v = atom->v, ** f = atom->f;
  int dim = domain->dimension;
  int m = 0;
  int last = first + n;
  for (int i=first; i<last; i++) {
    for (int k=0; k<dim; k++) x[i][k] = buf[m++];
    for (int k=0; k<dim; k++) v[i][k] = buf[m++];
    for (int k=0; k<dim; k++) f[i][k] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */
