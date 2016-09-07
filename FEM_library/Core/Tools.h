/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>

#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

/*************************************************************************/

double Mdot(int dim, double X[dim][dim], double Y[dim][dim]);


double ScalarProduct(int dim, double x[dim], double y[dim]);


void MatrixVector(int dim1, int dim2, double A[dim1][dim2], double x[dim2], double y[dim1]);


void MatrixSum(int dim, double X[dim][dim], double Y[dim][dim] );


void MatrixSumAlpha(int dim, double alpha, double X[dim][dim], double beta, double Y[dim][dim], double result[dim][dim] );


double MatrixDeterminant2(int dim, double A[dim][dim]);


double MatrixDeterminant3(int dim, double A[dim][dim]);


double MatrixDeterminant(int dim, double A[dim][dim]);


void MatrixInvT(int dim, double A[dim][dim], double invAT[dim][dim] );


void MatrixInvT3(int dim, double A[dim][dim], double invAT[dim][dim] );


void MatrixInv3(int dim, double A[dim][dim], double invA[dim][dim] );


void MatrixProduct(int dim, double X[dim][dim], double Y[dim][dim], double result[dim][dim] );


void MatrixScalar(int dim, double alpha, double X[dim][dim], double result[dim][dim] );


void MatrixProductAlpha(int dim, double alpha, double X[dim][dim], double Y[dim][dim], double result[dim][dim] );


void MatrixProductAlphaT1(int dim, double alpha, double X[dim][dim], double Y[dim][dim], double result[dim][dim] );


void MatrixProductAlphaT2(int dim, double alpha, double X[dim][dim], double Y[dim][dim], double result[dim][dim] );


void MatrixProductAlphaT3(int dim, double alpha, double X[dim][dim], double Y[dim][dim], double result[dim][dim] );


void MatrixProductQ1(int dim, int numQuadPoints, double X[dim][dim][numQuadPoints], double Y[dim][dim], double result[dim][dim], int q );


double Trace(int dim, double X[dim][dim]);


double TraceQ(int dim, int numQuadPoints, double X[dim][dim][numQuadPoints], int q);


void compute_GreenStrainTensor(int dim, int numQuadPoints, double F[dim][dim][numQuadPoints], double Id[dim][dim], double E[dim][dim][numQuadPoints], int q );


void compute_DerGreenStrainTensor(int dim, int numQuadPoints, double F[dim][dim][numQuadPoints], double dF[dim][dim], double dE[dim][dim], int q );


#endif