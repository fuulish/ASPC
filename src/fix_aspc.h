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

#ifdef FIX_CLASS

FixStyle(aspc,FixASPC)    // This registers this fix class with LAMMPS.

#else
*/

#ifndef LMP_FIX_ASPC_H
#define LMP_FIX_ASPC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixASPC : public Fix {
 public:
  FixASPC(class LAMMPS *, int, char **);
  ~FixASPC();
  void init();
  void reset_vectors();
  void sort_vectors();
  void update_history();
  // void allocate();
  // void deallocate();
  // void setup_style();
  void final_integrate();
  void initial_integrate(int);
  void pre_force(int);
  virtual void setup_pre_force(int);
  int modify_param(int narg, char **arg);

  // void setup_pre_force(int);
  // void setup_pre_force_respa(int, int);
  // void pre_force_respa(int, int, int);
  // void pre_force(int);

  // derived child classes must provide these functions
  // FUDO| we'll do it differently, we'll assume 

  // virtual int iterate(int) = 0;
  // virtual void setup_style() = 0;

  int setmask();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  void add_vector(int);
  double *request_vector(int);

 protected:
  int length, nord;             // length of the predictor
  int dim;                     // dimensionality of the thing to be aspc'ed
  int zerovels;                // zero-out particle velocities if needed
  double damp;                 // damping factor to be used in simulation, if set to something negative, default is used
  int what;                  // what to predict/correct
  int first;                   // at which timestep is the fix initialized
  bigint ndoftotal;           // total dof for entire problem

  enum{NONE,COORDS,QEQ,QEQREAX}; // possible things to predict/correct

  int nvec;                   // local atomic dof = length of xvec
  double *qty;                // variables for atomic dof, as 1d vector

  class Compute *pe_compute;        // compute for potential energy
  double ecurrent;                  // current potential energy
  double dmax;                // max dist to move any atom in one step

  //FU| copying in bulk - still needs sorting

  int triclinic;              // 0 if domain is orthog, 1 if triclinic
  int eflag,vflag;            // flags for energy/virial computation
  int virial_style;           // compute virial explicitly or implicitly
  int external_force_clear;   // clear forces locally or externally

  int torqueflag,extraflag;

  int pair_compute_flag;            // 0 if pair->compute is skipped
  int kspace_compute_flag;          // 0 if kspace->compute is skipped

  int neigh_every,neigh_delay,neigh_dist_check;  // neighboring params

  void force_clear();

  // void setup();
  // void setup_minimal();

  void generate_coefficients(int);
  double get_Bj(int n, int k);
  virtual void predict();
  virtual void correct();
  virtual void cpy2hist();
  void reset_history();
  void clear_non_group();

  //FUDO| double-check this function
  long long combi(int n,int k);

  int tstart;
  int nonlyhist;
  int filledhist;
  double *coeffs;

 private:
  int nvector;
  int *peratom;
  double **vectors;
  double bzr;
};

}

#endif
// #endif
