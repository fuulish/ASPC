/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* FUDO| this should not be done, because we only want to use the derived classes
*/

#ifdef FIX_CLASS

FixStyle(aspc/drude,FixASPCDrude)    // This registers this fix class with LAMMPS.

#else

#ifndef LMP_FIX_ASPC_DRUDE_H
#define LMP_FIX_ASPC_DRUDE_H

#include "fix.h"
#include "fix_drude.h"
#include "fix_aspc.h"

namespace LAMMPS_NS {

class FixASPCDrude : public FixASPC {
 public:
  FixASPCDrude(class LAMMPS *, int, char **);
  ~FixASPCDrude();
  void init();
  int modify_param(int narg, char **arg);
  int check_convergence();
  double calc_spring_forces_energy();
  void predict();
  void cpy2hist();
  void setup_pre_force(int vflag);
  void r2r_indices_forward();
  void r2r_indices_reverse();

 protected:
  void correct();
  FixDrude * fix_drude;
  tagint *drudeid;
  double ftol;
  int scf, nfail, faild, neval, printconv;
  double *dx, *hrm;

 private:
  double kd;
  char *id_ef;
  class Compute *c_ef;
};

}

#endif
#endif
