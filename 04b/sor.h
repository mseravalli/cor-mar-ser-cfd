#ifndef __SOR_H_
#define __SOR_H_

/**
 * One GS iteration for the pressure Poisson equation. Besides, the routine must 
 * also set the boundary values for P according to the specification. The 
 * residual for the termination criteria has to be stored in res.
 * 
 * An \omega = 1 GS - implementation is given within sor.c.
 */
void sor(
  double omg,
  double dx,
  double dy,
  int imax,
  int jmax,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  int chunk,
  double* bufSend,
  double* bufRecv,
  double **P,
  double **RS,
  double *res
);


#endif
