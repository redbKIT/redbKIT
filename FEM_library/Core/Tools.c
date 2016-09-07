/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "Tools.h"

/*************************************************************************/

double Mdot(int dim, double X[dim][dim], double Y[dim][dim])
{
    int d1, d2;
    double Z = 0;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Z = Z + X[d1][d2] * Y[d1][d2];
        }
    }
    return Z;
}
/*************************************************************************/
double ScalarProduct(int dim, double x[dim], double y[dim])
{
    double result = 0;
    int d1;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        result += x[d1] * y[d1];
    }
    return result;
}
/*************************************************************************/
void MatrixVector(int dim1, int dim2, double A[dim1][dim2], double x[dim2], double y[dim1])
{
    int d1, d2;
    for (d1 = 0; d1 < dim1; d1 = d1 + 1 )
    {
        y[d1] = 0;
        for (d2 = 0; d2 < dim2; d2 = d2 + 1 )
        {
            y[d1] += A[d1][d2] * x[d2];
        }
    }    
}

/*************************************************************************/
void MatrixSum(int dim, double X[dim][dim], double Y[dim][dim] )
{
    int d1, d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            X[d1][d2] = X[d1][d2] + Y[d1][d2];
        }
    }
}
/*************************************************************************/
void MatrixSumAlpha(int dim, double alpha, double X[dim][dim], double beta, double Y[dim][dim], double result[dim][dim] )
{
    int d1, d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            result[d1][d2] = alpha * X[d1][d2] + beta * Y[d1][d2];
        }
    }
}
/*************************************************************************/
double MatrixDeterminant2(int dim, double A[dim][dim])
{
    return A[0][0]*A[1][1]-A[0][1]*A[1][0];
}
/*************************************************************************/
double MatrixDeterminant3(int dim, double A[dim][dim])
{
    return A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2]) -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])  +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
}
/*************************************************************************/
double MatrixDeterminant(int dim, double A[dim][dim])
{
    double determinant = 0;
    if ( dim == 2 )
    {
        determinant = MatrixDeterminant2(2, A);
    }
    
    if ( dim == 3 )
    {
        determinant = MatrixDeterminant3(3, A);
    }
        
    return determinant;
}
/*************************************************************************/
void MatrixInvT(int dim, double A[dim][dim], double invAT[dim][dim] )
{
    if ( dim == 2 )
    {
        double det = MatrixDeterminant2(2, A);
        
        double invdet = 1/det;
        
        invAT[0][0] =   A[1][1]*invdet;
        invAT[1][0] =  -A[0][1]*invdet;
        invAT[0][1] =  -A[1][0]*invdet;
        invAT[1][1] =   A[0][0]*invdet;
    }
    
    if ( dim == 3 )
    {
        double det = MatrixDeterminant3(3, A);
        
        double invdet = 1/det;
        
        invAT[0][0] =  (A[1][1]*A[2][2]-A[2][1]*A[1][2])*invdet;
        invAT[1][0] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
        invAT[2][0] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
        invAT[0][1] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
        invAT[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
        invAT[2][1] = -(A[0][0]*A[1][2]-A[1][0]*A[0][2])*invdet;
        invAT[0][2] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
        invAT[1][2] = -(A[0][0]*A[2][1]-A[2][0]*A[0][1])*invdet;
        invAT[2][2] =  (A[0][0]*A[1][1]-A[1][0]*A[0][1])*invdet;
    }
    
}
/*************************************************************************/
void MatrixInvT3(int dim, double A[dim][dim], double invAT[dim][dim] )
{
   
    double det = MatrixDeterminant3(3, A);
    
    double invdet = 1/det;
    
    invAT[0][0] =  (A[1][1]*A[2][2]-A[2][1]*A[1][2])*invdet;
    invAT[1][0] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
    invAT[2][0] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
    invAT[0][1] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
    invAT[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
    invAT[2][1] = -(A[0][0]*A[1][2]-A[1][0]*A[0][2])*invdet;
    invAT[0][2] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
    invAT[1][2] = -(A[0][0]*A[2][1]-A[2][0]*A[0][1])*invdet;
    invAT[2][2] =  (A[0][0]*A[1][1]-A[1][0]*A[0][1])*invdet;

}
/*************************************************************************/
void MatrixInv3(int dim, double A[dim][dim], double invA[dim][dim] )
{
   
    double det = MatrixDeterminant3(3, A);
    
    double invdet = 1/det;
    
    invA[0][0] =  (A[1][1]*A[2][2]-A[2][1]*A[1][2])*invdet;
    invA[0][1] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
    invA[0][2] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
    invA[1][0] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
    invA[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
    invA[1][2] = -(A[0][0]*A[1][2]-A[1][0]*A[0][2])*invdet;
    invA[2][0] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
    invA[2][1] = -(A[0][0]*A[2][1]-A[2][0]*A[0][1])*invdet;
    invA[2][2] =  (A[0][0]*A[1][1]-A[1][0]*A[0][1])*invdet;

}
/*************************************************************************/
void MatrixProduct(int dim, double X[dim][dim], double Y[dim][dim], double result[dim][dim] )
{
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            result[d1][d2] = 0;
            for (d3 = 0; d3 < dim; d3 = d3 + 1 )
            {
                result[d1][d2] = result[d1][d2] + X[d1][d3]*Y[d3][d2];
            }
        }
    }
}
/*************************************************************************/
void MatrixScalar(int dim, double alpha, double X[dim][dim], double result[dim][dim] )
{
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            result[d1][d2] = alpha * X[d1][d2];
        }
    }
}
/*************************************************************************/
void MatrixProductAlpha(int dim, double alpha, double X[dim][dim], double Y[dim][dim], double result[dim][dim] )
{
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            result[d1][d2] = 0;
            for (d3 = 0; d3 < dim; d3 = d3 + 1 )
            {
                result[d1][d2] = result[d1][d2] + X[d1][d3]*Y[d3][d2];
            }
            result[d1][d2] = alpha * result[d1][d2];
        }
    }
}
/*************************************************************************/
void MatrixProductAlphaT1(int dim, double alpha, double X[dim][dim], double Y[dim][dim], double result[dim][dim] )
{
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            result[d1][d2] = 0;
            for (d3 = 0; d3 < dim; d3 = d3 + 1 )
            {
                result[d1][d2] = result[d1][d2] + X[d3][d1]*Y[d3][d2];
            }
            result[d1][d2] = alpha * result[d1][d2];
        }
    }
}
/*************************************************************************/
void MatrixProductAlphaT2(int dim, double alpha, double X[dim][dim], double Y[dim][dim], double result[dim][dim] )
{
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            result[d1][d2] = 0;
            for (d3 = 0; d3 < dim; d3 = d3 + 1 )
            {
                result[d1][d2] = result[d1][d2] + X[d1][d3]*Y[d2][d3];
            }
            result[d1][d2] = alpha * result[d1][d2];
        }
    }
}
/*************************************************************************/
void MatrixProductAlphaT3(int dim, double alpha, double X[dim][dim], double Y[dim][dim], double result[dim][dim] )
{
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            result[d1][d2] = 0;
            for (d3 = 0; d3 < dim; d3 = d3 + 1 )
            {
                result[d1][d2] = result[d1][d2] + X[d3][d1]*Y[d2][d3];
            }
            result[d1][d2] = alpha * result[d1][d2];
        }
    }
}
/*************************************************************************/
void MatrixProductQ1(int dim, int numQuadPoints, double X[dim][dim][numQuadPoints], double Y[dim][dim], double result[dim][dim], int q )
{
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            result[d1][d2] = 0;
            for (d3 = 0; d3 < dim; d3 = d3 + 1 )
            {
                result[d1][d2] = result[d1][d2] + X[d1][d3][q]*Y[d3][d2];
            }
        }
    }
}
/*************************************************************************/
double Trace(int dim, double X[dim][dim])
{
    double T = 0;
    int d1;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        T = T + X[d1][d1];
    }
    return T;
}
/*************************************************************************/
double TraceQ(int dim, int numQuadPoints, double X[dim][dim][numQuadPoints], int q)
{
    double T = 0;
    int d1;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        T = T + X[d1][d1][q];
    }
    return T;
}
/*************************************************************************/
void compute_GreenStrainTensor(int dim, int numQuadPoints, double F[dim][dim][numQuadPoints], double Id[dim][dim], double E[dim][dim][numQuadPoints], int q )
{
    
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            double tmp = 0;
            for (d3 = 0; d3 < dim; d3 = d3 + 1 )
            {
                tmp = tmp + F[d3][d1][q] * F[d3][d2][q];
            }
            E[d1][d2][q] = 0.5 * ( tmp - Id[d1][d2] );
        }
    }
}
/*************************************************************************/
void compute_DerGreenStrainTensor(int dim, int numQuadPoints, double F[dim][dim][numQuadPoints], double dF[dim][dim], double dE[dim][dim], int q )
{
    
    int d1, d2, d3;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            double tmp1 = 0;
            double tmp2 = 0;
            for (d3 = 0; d3 < dim; d3 = d3 + 1 )
            {
                tmp1 = tmp1 + dF[d3][d1] * F[d3][d2][q];
                tmp2 = tmp2 + F[d3][d1][q]  * dF[d3][d2];
            }
            dE[d1][d2] = 0.5 * ( tmp1 + tmp2 );
        }
    }
}
/*************************************************************************/
