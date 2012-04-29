#ifndef __DISC_H__                                                                                    
#define __DISC_H__


/******** DERIVATIVES FOR EQUATION (1) ********/
double** du2dx(double** U, double dx, int nrl, int nrh, int ncl, int nch);
double** duvdy(double** U, double **V, double dy, int nrl, int nrh, int ncl, int nch);
double** d2udx2(double** U, double dx, int nrl, int nrh, int ncl, int nch); 
double** d2udy2(double** U, double dy, int nrl, int nrh, int ncl, int nch); 
double** dpdx(double** P, double dx, int nrl, int nrh, int ncl, int nch);


/******** DERIVATIVES FOR EQUATION (2) ********/
double** dv2dy(double** V, double dy, int nrl, int nrh, int ncl, int nch);
double** duvdx(double** U, double **V, double dx, int nrl, int nrh, int ncl, int nch);
double** d2vdx2(double** V, double dx, int nrl, int nrh, int ncl, int nch); 
double** d2vdy2(double** V, double dy, int nrl, int nrh, int ncl, int nch); 
double** dpdy(double** P, double dy, int nrl, int nrh, int ncl, int nch);

#endif




