#ifndef __MATRIX_OP_H__                                                                                  
#define __MATRIX_OP_H__ 

/*adds scalar to all elements of a matrix*/
double** add_scalar(double** M, int nrl, int nrh, int ncl, int nch, double s);

/*multiplies all the elements of a matrix by a scalar*/
double** mult_scalar(double** M, int nrl, int nrh, int ncl, int nch, double s);

/*adds two matrices with same dimension */
double** add_mat(double** M, double** N, int nrl, int nrh, int ncl, int nch);

/*subtracts the second matrix from the first one with same dimension */
double** sub_mat(double** M, double** N, int nrl, int nrh, int ncl, int nch);

#endif
