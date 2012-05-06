#include "matrix_op.h"
#include "helper.h"


double** add_scalar(double** M, double s, int nrl, int nrh, int ncl, int nch, double** result)
{
    int i = 0;
    int j = 0;
    
    for(i = nrl; i <= nrh; ++i)
    {
        for(j = ncl; j <= nch; ++j)
        {
            result[i][j] = M[i][j] + s;
        }
    }
    
    return result;
}

double** mult_scalar(double** M, double s, int nrl, int nrh, int ncl, int nch, double** result)
{
    int i = 0;
    int j = 0;
    
    for(i = nrl; i <= nrh; ++i)
    {
        for(j = ncl; j <= nch; ++j)
        {
            result[i][j] = M[i][j] * s;
        }
    }
    
    return result;
}


double** add_mat(double** M, double** N, int nrl, int nrh, int ncl, int nch, double** result)
{
    int i = 0;
    int j = 0;
    
    for(i = nrl; i <= nrh; ++i)
    {
        for(j = ncl; j <= nch; ++j)
        {
            result[i][j] = M[i][j] + N[i][j];
        }
    }
       
    return result;
}


double** sub_mat(double** M, double** N, int nrl, int nrh, int ncl, int nch, double** result)
{
    int i = 0;
    int j = 0;
    
    for(i = nrl; i <= nrh; ++i)
    {
        for(j = ncl; j <= nch; ++j)
        {
            result[i][j] = M[i][j] - N[i][j];
        }
    }
       
    return result;
}
