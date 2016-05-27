/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]
#ifdef _OPENMP
#include <omp.h>
#else
#warning "OpenMP not enabled. Compile with mex CSM_ComputeStress_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

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

void StVenantKirchhoffMaterial(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(noe,dim*dim, mxREAL);
    
    double* Sigma    = mxGetPr(plhs[0]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);
    
    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);
    
    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[4]);
    
    double GradUh[dim][dim][NumQuadPoints];
    
    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }
    
    double F[dim][dim][NumQuadPoints];
    double E[dim][dim][NumQuadPoints];
    double P_Uh[dim][dim];
        
    double* material_param = mxGetPr(prhs[2]);
    double Young = material_param[0];
    double Poisson = material_param[1];
    double mu = Young / (2 + 2 * Poisson);
    double lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,Sigma,U_h) private(gradphi,F,E,P_Uh,GradUh,ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        double traceE[NumQuadPoints];
        q = 0;
        
        /* Compute Gradient of Basis functions*/
        for (k = 0; k < nln; k = k + 1 )
        {
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                gradphi[d1][k][q] = 0;
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    gradphi[d1][k][q] = gradphi[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                }
            }
        }
        
        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                GradUh[d1][d2][q] = 0;
                for (k = 0; k < nln; k = k + 1 )
                {
                    int e_k;
                    e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                    GradUh[d1][d2][q] = GradUh[d1][d2][q] + U_h[e_k] * gradphi[d2][k][q];
                }
                F[d1][d2][q]  = Id[d1][d2] + GradUh[d1][d2][q];
            }
        }
        compute_GreenStrainTensor(dim, NumQuadPoints, F, Id, E, q );
        traceE[q] = TraceQ(dim, NumQuadPoints, E, q);
        
        double P1[dim][dim];
        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                P1[d1][d2] =  ( 2 * mu * E[d1][d2][q] + lambda * traceE[q] * Id[d1][d2] );
            }
        }
        MatrixProductQ1(dim, NumQuadPoints, F, P1, P_Uh, q);
        
        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                /*Sigma[ie][d1*dim+d2] =  P_Uh[d1][d2] ;*/
                Sigma[ie+(d1+d2*dim)*noe] =  P_Uh[d1][d2] ;
            }
        }
    }
    
}
/*************************************************************************/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=11) {
        mexErrMsgTxt("11 inputs are required.");
    } else if(nlhs>1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    char *Material_Model = mxArrayToString(prhs[1]);
    
    /*
    if (strcmp(Material_Model, "Linear")==0)
    {
            LinearElasticMaterial(plhs, prhs);
    }
    */

    if (strcmp(Material_Model, "StVenantKirchhoff")==0)
    {
            StVenantKirchhoffMaterial(plhs, prhs);
    }
    
    mxFree(Material_Model);
}
/*************************************************************************/

