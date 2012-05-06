#ifndef __MATRIX_OP_H__                                                                                  
#define __MATRIX_OP_H__ 

/*adds scalar to all elements of a matrix*/
void add_scalar(double** M, double s, int nrl, int nrh, int ncl, int nch, double** result);

/*multiplies all the elements of a matrix by a scalar*/
void mult_scalar(double** M, double s, int nrl, int nrh, int ncl, int nch, double** result);

/*adds two matrices with same dimension */
void add_mat(double** M, double** N, int nrl, int nrh, int ncl, int nch, double** result);

/*subtracts the second matrix from the first one with same dimension */
void sub_mat(double** M, double** N, int nrl, int nrh, int ncl, int nch, double** result);

#endif
