/* ----------------------------------------------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------------------------
   Contributing Author: Utkarsh Kapoor (Texas A&M), Young Chan Kim (NRL Navy), Jeetain Mittal (Texas A&M)
------------------------------------------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj1210dna_2bead.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLj1210DNA_2Bead::PairLj1210DNA_2Bead(LAMMPS *lmp) : Pair(lmp) 
{
  TWO_1_3 = pow(2.0,(1.0/3.0));
}

/* ---------------------------------------------------------------------- */

PairLj1210DNA_2Bead::~PairLj1210DNA_2Bead()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
 
    memory->destroy(a);
    memory->destroy(b);
    memory->destroy(c);
    memory->destroy(d);
    memory->destroy(epsilon); 
    memory->destroy(sigma);
    memory->destroy(q);
    memory->destroy(qc);
    memory->destroy(offset);

    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(lj5);
    memory->destroy(lj6);
    memory->destroy(lj7);
    memory->destroy(lj8);
    memory->destroy(lj9);
    memory->destroy(lj10);
    memory->destroy(lj11);
    memory->destroy(lj12);
  }
}

/* ---------------------------------------------------------------------- */

void PairLj1210DNA_2Bead::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,k,m;
  int imol, jmol, itag, jtag,ktag,mtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double delr1[3], delr2[3], r1, r2, rsq1, rsq2,rsq,r2inv,r6inv,r4inv,forcesthb,factor_lj;
  double forcelj;
  double ac, PI, angle;
  PI = 3.1415926535;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *molecule = atom->molecule;           //point to molecule ID
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {


    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];       
    itag = tag[i];       
    imol = molecule[i];           //get the molecule id for i 
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      jtag = tag[j];
      jmol = molecule[j];         //get the molecule id for j 

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r4inv = r2inv*r2inv;
        r6inv = r2inv*r2inv*r2inv;

        // if in the same strand && adjcent to each other (1-2 interactions)
        if ((imol == jmol) && ((jtag == itag + 2) || (jtag == itag - 2))) {
                 if (rsq < qc[itype][jtype]*qc[itype][jtype]) {
                      forcesthb = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]*r4inv);
                    }
                 else forcesthb = 0;
         }

        // if in the same strand && 1-3 interactions
        else if ((imol == jmol) && ((jtag == itag + 4) || (jtag == itag - 4))) {
                 if (rsq < TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) {
                      forcesthb = r6inv*(lj5[itype][jtype]*r6inv-lj6[itype][jtype]);
                 } 
                 else forcesthb = 0.0;
        }

        // 1-4, 1-5 ... and intermolecular interations
        else if (q[itype][jtype] == 1) { // if stacking is on
                 if (rsq < TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) {
                    forcesthb = r6inv*(lj5[itype][jtype]*r6inv-lj6[itype][jtype]);
                 }
                 else forcesthb = 0.0;
        }

        else { // if HB is on
               forcesthb = r6inv * (lj9[itype][jtype]*r6inv - lj10[itype][jtype]*r4inv);
        }
        
        fpair = factor_lj*forcesthb*r2inv;
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          // if in the same strand && adjcent to each other (1-2 interactions)
          if ((imol == jmol) && ((jtag == itag + 2) || (jtag == itag - 2))) {        
                  if (rsq < qc[itype][jtype]*qc[itype][jtype]) {
                      evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]*r4inv);
                  }
                  else evdwl = 0;	
          }

          // if in the same strand && 1-3 interactions
          else if ((imol == jmol) && ((jtag == itag + 4) || (jtag == itag - 4))) {
                  if (rsq < TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) {
                      evdwl = r6inv*(lj7[itype][jtype]*r6inv-lj8[itype][jtype]) - offset[itype][jtype];
                  }
                  else evdwl =0;
          }

          // 1-4, 1-5 ... and intermolecular interations
          else if (q[itype][jtype] == 1) { // if stacking is on
                  if (rsq < TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) {
                      evdwl = r6inv*(lj7[itype][jtype]*r6inv-lj8[itype][jtype]) - offset[itype][jtype];
                  }
                  else evdwl =0;
          }

          else { // if HB is on
                  evdwl = (r6inv*(lj11[itype][jtype]*r6inv-lj12[itype][jtype]*r4inv));
          }
    
          evdwl *= factor_lj;
        }
      if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLj1210DNA_2Bead::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");

  memory->create(a,n+1,n+1,"pair:a");
  memory->create(b,n+1,n+1,"pair:b");
  memory->create(c,n+1,n+1,"pair:c");
  memory->create(d,n+1,n+1,"pair:d");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(q,n+1,n+1,"pair:q");
  memory->create(qc,n+1,n+1,"pair:qc");
  memory->create(offset,n+1,n+1,"pair:offset");

  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(lj5,n+1,n+1,"pair:lj5");
  memory->create(lj6,n+1,n+1,"pair:lj6");
  memory->create(lj7,n+1,n+1,"pair:lj7");
  memory->create(lj8,n+1,n+1,"pair:lj8");
  memory->create(lj9,n+1,n+1,"pair:lj9");
  memory->create(lj10,n+1,n+1,"pair:lj10");
  memory->create(lj11,n+1,n+1,"pair:lj11");
  memory->create(lj12,n+1,n+1,"pair:lj12");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLj1210DNA_2Bead::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");
  cut_global = utils::numeric(FLERR,arg[0],false,lmp); //UK

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLj1210DNA_2Bead::coeff(int narg, char **arg)
{
  if (narg < 10 || narg > 11) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);//UK
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);//UK  

  double a_one = utils::numeric(FLERR,arg[2],false,lmp);//UK
  double b_one = utils::numeric(FLERR,arg[3],false,lmp);//UK
  double c_one = utils::numeric(FLERR,arg[4],false,lmp);//UK
  double d_one = utils::numeric(FLERR,arg[5],false,lmp);//UK
  double epsilon_one = utils::numeric(FLERR,arg[6],false,lmp);//UK
  double sigma_one = utils::numeric(FLERR,arg[7],false,lmp);//UK
  double q_one = utils::numeric(FLERR,arg[8],false,lmp);//UK
  double qc_one = utils::numeric(FLERR,arg[9],false,lmp);//UK
  double cut_one = cut_global;
  if (narg == 11) cut_one = utils::numeric(FLERR,arg[10],false,lmp); //UK
 
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a[i][j] = a_one;
      b[i][j] = b_one;
      c[i][j] = c_one;
      d[i][j] = d_one;
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      q[i][j] = q_one;
      qc[i][j] = qc_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLj1210DNA_2Bead::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  lj1[i][j] = 60.0 * a[i][j] * pow(b[i][j], 12.0);
  lj2[i][j] = 60.0 * a[i][j] * pow(b[i][j], 10.0);
  lj3[i][j] = 5.0 * a[i][j] * pow(b[i][j], 12.0);
  lj4[i][j] = 6.0 * a[i][j] * pow(b[i][j], 10.0);
  lj5[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj6[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj7[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj8[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj9[i][j] = 60.0 * c[i][j] * pow(d[i][j], 12.0);
  lj10[i][j] = 60.0 * c[i][j] * pow(d[i][j], 10.0);
  lj11[i][j] = 5.0 * c[i][j] * pow(d[i][j], 12.0);
  lj12[i][j] = 6.0 * c[i][j] * pow(d[i][j], 10.0);

  
  offset[i][j] = -epsilon[i][j];
 
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  lj5[j][i] = lj5[i][j];
  lj6[j][i] = lj6[i][j];
  lj7[j][i] = lj7[i][j];
  lj8[j][i] = lj8[i][j];
  lj9[j][i] = lj9[i][j];
  lj10[j][i] = lj10[i][j];
  lj11[j][i] = lj11[i][j];
  lj12[j][i] = lj12[i][j];
     

  a[j][i] = a[i][j];
  b[j][i] = b[i][j];
  c[j][i] = c[i][j];
  d[j][i] = d[i][j];
  q[j][i] = q[i][j];
  qc[j][i] = qc[i][j];
  sigma[j][i] = sigma[i][j];
  epsilon[j][i] = epsilon[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction 
  // count total # of atoms of type I and J via Allreduce 

  if (tail_flag) { 
     int *type = atom->type; 
     int nlocal = atom->nlocal; 

     double count[2],all[2]; 
     count[0] = count[1] = 0.0; 
     for (int k = 0; k < nlocal; k++) { 
       if (type[k] == i) count[0] += 1.0; 
       if (type[k] == j) count[1] += 1.0; 
     } 
     MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

   } 

   return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLj1210DNA_2Bead::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&a[i][j],sizeof(double),1,fp);
        fwrite(&b[i][j],sizeof(double),1,fp);
	fwrite(&c[i][j],sizeof(double),1,fp);
	fwrite(&d[i][j],sizeof(double),1,fp);
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
	fwrite(&q[i][j],sizeof(double),1,fp);
	fwrite(&qc[i][j],sizeof(double),1,fp);
	fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLj1210DNA_2Bead::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
          utils::sfread(FLERR,&a[i][j],sizeof(double),1,fp,NULL,error); //UK
          utils::sfread(FLERR,&b[i][j],sizeof(double),1,fp,NULL,error); //UK
          utils::sfread(FLERR,&c[i][j],sizeof(double),1,fp,NULL,error); //UK
          utils::sfread(FLERR,&d[i][j],sizeof(double),1,fp,NULL,error); //UK
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,NULL,error); //UK
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,NULL,error); //UK
          utils::sfread(FLERR,&q[i][j],sizeof(double),1,fp,NULL,error); //UK
          utils::sfread(FLERR,&qc[i][j],sizeof(double),1,fp,NULL,error); //UK
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,NULL,error); //UK
         
	}
	MPI_Bcast(&a[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&b[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&c[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&d[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&q[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&qc[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLj1210DNA_2Bead::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLj1210DNA_2Bead::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,NULL,error); //UK
    utils::sfread(FLERR,&offset_flag,sizeof(double),1,fp,NULL,error); //UK
    utils::sfread(FLERR,&mix_flag,sizeof(double),1,fp,NULL,error); //UK
    utils::sfread(FLERR,&tail_flag,sizeof(double),1,fp,NULL,error); //UK
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairLj1210DNA_2Bead::single(int i, int j, int itype, int jtype,
			double rsq, double factor_coul, double factor_lj,
			double &fforce)
{
  double r2inv,r6inv,r4inv,forcesthb,phisthb,forcelj;
  int imol, jmol, itag, jtag;

  imol = atom->molecule[i];
  jmol = atom->molecule[j];
  itag = atom->tag[i];
  jtag = atom->tag[j];

  r2inv = 1.0/rsq;
  r4inv = r2inv*r2inv;
  r6inv = r2inv*r2inv*r2inv;

          // if in the same strand && adjcent to each other (1-2 interactions)
          if ((imol == jmol) && ((jtag == itag + 2) || (jtag == itag - 2))) {
		  if (rsq < qc[itype][jtype]*qc[itype][jtype]) {
  	              forcesthb = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]*r4inv);
		  }
          }

          // if in the same strand && 1-3 interactions
          else if ((imol == jmol) && ((jtag == itag + 4) || (jtag == itag - 4))) {
		if (rsq < TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) {
		          forcesthb = r6inv*(lj5[itype][jtype]*r6inv-lj6[itype][jtype]);
		 }
	  }

          // 1-4, 1-5 ... and intermolecular interations
          else if (q[itype][jtype] == 1) { // if stacking is on
		  if (rsq < TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) {
			  forcesthb = r6inv*(lj5[itype][jtype]*r6inv-lj6[itype][jtype]);
		  }
	  }

	  else { // if HB is on
                 forcesthb = (r6inv * (lj9[itype][jtype]*r6inv - lj10[itype][jtype]*r4inv));
          }          
        
        
  fforce = factor_lj*forcesthb*r2inv;

 
          // if in the same strand && adjcent to each other (1-2 interactions)
            if ((imol == jmol) && ((jtag == itag + 2) || (jtag == itag - 2))) {
		     if (rsq < qc[itype][jtype]*qc[itype][jtype]) {
	                phisthb = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]*r4inv);
		     }
            }

          // if in the same strand && 1-3 interactions
            else if ((imol == jmol) && ((jtag == itag + 4) || (jtag == itag - 4))) {
		    if (rsq < TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) {
	                phisthb = r6inv*(lj7[itype][jtype]*r6inv-lj8[itype][jtype]) - offset[itype][jtype];
		    }
            }

          // 1-4, 1-5 ... and intermolecular interations
            else if (q[itype][jtype] == 1) { // if stacking is on
		    if (rsq < TWO_1_3*sigma[itype][jtype]*sigma[itype][jtype]) {
			    phisthb = r6inv*(lj7[itype][jtype]*r6inv-lj8[itype][jtype]) - offset[itype][jtype];
		    }
            }

            else { // if HB is on
                    phisthb = (r6inv*(lj11[itype][jtype]*r6inv-lj12[itype][jtype]*r4inv));
            }
        
  return factor_lj*phisthb;
}

/* ---------------------------------------------------------------------- */


