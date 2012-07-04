#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  int kmax,
  double dx,
  double dy,
  int wl,
  int wr,
  int wt,
  int wb,
  double** U,
  double** V,
  double** F,
  double** G,
  double** P,
  double***C,
  int **Flag,
  int **Sources,
  double* K
);

/*Set special boundary*/

void spec_boundary_val(
    char *problem,
    int imax,
    int jmax,
    double dx,
    double dy,
    double Re,
    double deltaP,
    double **U,
    double **V,
    double **P,
    double*** C

);

#endif
