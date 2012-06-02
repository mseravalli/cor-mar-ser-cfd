#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int il,
  int ir,
  int jt,
  int jb,
  int omg_i,
  int omg_j,
  int iproc,
  int jproc,
  double **U,
  double **V
);

#endif
