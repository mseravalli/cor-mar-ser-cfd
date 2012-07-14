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
  int wb,
  int wt,
  int cl,
  int cr,
  int cb,
  int ct,
  double** U,
  double** V,
  double** F,
  double** G,
  double** P,
  double***C,
  double** T,
  int **Flag,
  int **Sources,
  double* C0
);

/*Set special boundary*/

void spec_boundary_val(
    char *problem,
    int imax,
    int jmax,
    int kmax,
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
