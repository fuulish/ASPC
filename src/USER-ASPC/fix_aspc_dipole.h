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

FixStyle(aspc/dipole,FixASPCDipole)    // This registers this fix class with LAMMPS.

#else

#ifndef LMP_FIX_ASPC_DIPOLE_H
#define LMP_FIX_ASPC_DIPOLE_H

#include "fix.h"
#include "fix_aspc.h"

namespace LAMMPS_NS {

class FixASPCDipole : public FixASPC {
 public:
  FixASPCDipole(class LAMMPS *, int, char **);
  ~FixASPCDipole();
  void init();
  int modify_param(int narg, char **arg);
  int check_convergence(double, double, double **, double *);
  double compute_pe_scalar();

 protected:
  void correct();
  void predict();
  double energy_force_es(int);
  double ftol;
  int scf, neval, printconv;
  double *dx, *oldf;
  int ndx_dim;                     // actual dimensionality used for indexing the thing to be aspc'ed

 private:
  double fpieps;
};

}

#endif
#endif
