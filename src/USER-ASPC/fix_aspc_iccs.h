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

FixStyle(aspc/iccs,FixASPCICCS)    // This registers this fix class with LAMMPS.

#else

#ifndef LMP_FIX_ASPC_ICCS_H
#define LMP_FIX_ASPC_ICCS_H

#include "fix.h"
#include "fix_aspc.h"

namespace LAMMPS_NS {

class FixASPCICCS : public FixASPC {
 public:
  FixASPCICCS(class LAMMPS *, int, char **);
  ~FixASPCICCS();
  void init();
  int modify_param(int narg, char **arg);
  void predict();
  void setup_pre_force(int vflag);

  int pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc);
  void unpack_forward_comm(int n, int first, double *buf);
  void reset_vectors();

 protected:
  void correct();
  void iterate();
  void backup_charges();
  void calculate_charges_iccs();
  void initialize_charges();
  void update_charges();
  int check_convergence();
  void calculate_contrast();
  int scf, nfail, faild, neval, printconv;
  double *dx, *hrm;

 private:
  double conv;
  char *id_ef, *id_diel, *id_area, *id_srfx, *id_srfy, *id_srfz;
  class Compute *c_ef;
  int qinit;

  double bulk_perm;
  double *p_diel, *p_area, *p_srfx, *p_srfy, *p_srfz;
  double *contrast, *qprv, *qnxt;
};

}

#endif
#endif
