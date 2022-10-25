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

#ifdef PAIR_CLASS

PairStyle(lj1210dna/3bead,PairLj1210DNA_3Bead)

#else

#ifndef PAIR_LJ1210_DNA_3BEAD_H
#define PAIR_LJ1210_DNA_3BEAD_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLj1210DNA_3Bead : public Pair {
 public:
  PairLj1210DNA_3Bead(class LAMMPS *);
  virtual ~PairLj1210DNA_3Bead();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);


 protected:
  double cut_global;
  double **cut;
  double **a,**b,**sigma,**epsilon,**q;
  double **offset;
  double **lj1,**lj2,**lj3,**lj4,**lj5,**lj6,**lj7,**lj8;
  double TWO_1_3;
  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
