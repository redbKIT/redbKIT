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
    #warning "OpenMP not enabled. Compile with mex ADR_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=15) {
        mexErrMsgTxt("15 inputs are required.");
    } else if(nlhs>6) {
        mexErrMsgTxt("Too many output arguments.");
    }

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    /**/
    plhs[0] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
       
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    double* myMcoef    = mxGetPr(plhs[3]);
    double* myRrows    = mxGetPr(plhs[4]);
    double* myRcoef    = mxGetPr(plhs[5]);
    
    /* copy the string data from prhs[0] into a C string input_ buf.    */
    char *OP_string = mxArrayToString(prhs[1]);
    int OP[4] = {0, 0, 0, 0};
    if (strcmp(OP_string, "diffusion")==0)
    {
        OP[0] = 1;
    }
    
    if (strcmp(OP_string, "transport")==0)
    {
        OP[1] = 1;
    }
    
    if (strcmp(OP_string, "reaction")==0)
    {
        OP[2] = 1;
    }
    
    if (strcmp(OP_string, "source")==0)
    {
        OP[3] = 1;
    }
    
    if (strcmp(OP_string, "all")==0)
    {
        OP[0] = 1;
        OP[1] = 1;
        OP[2] = 1;
        OP[3] = 1;
    }
    mxFree(OP_string);
    
    double C_t[dim];
    double C_d[dim][dim];
    
    double* TC_d   = mxGetPr(prhs[2]);
    double* TC_t   = mxGetPr(prhs[3]);
    
    int k,l;

    for (k = 0; k < dim; k = k + 1 )
    {
        for (l = 0; l < dim; l = l + 1 )
        {
            C_d[k][l] = 0;
        }
        C_t[k] = 0;
    }
    
    if ((int)(TC_d[0])==10 && (int)(TC_d[1])==10)
    {
        for (l = 0; l < dim; l = l + 1 )
        {
            C_d[l][l] = 1;
        }
    }
    else
    {
        C_d[(int)(TC_d[0]-1)][(int)(TC_d[1]-1)] = 1;
    }
    
    if ((int)(TC_t[0])==10)
    {
        for (l = 0; l < dim; l = l + 1 )
        {
            C_t[l] = 1;
        }
    }
    else
    {
        C_t[(int)(TC_t[0]-1)] = 1;
    }
    
    /* Local mass matrix (computed only once) with quadrature nodes */
    double LocalMass[nln][nln];
    int q;
    int NumQuadPoints     = mxGetN(prhs[10]);
    
    double* mu   = mxGetPr(prhs[6]);
    double* conv_field   = mxGetPr(prhs[7]);
    double* si   = mxGetPr(prhs[8]);
    double* f    = mxGetPr(prhs[9]);
    double* w   = mxGetPr(prhs[10]);
    double* invjac = mxGetPr(prhs[11]);
    double* detjac = mxGetPr(prhs[12]);
    double* phi = mxGetPr(prhs[13]);
    double* gradrefphi = mxGetPr(prhs[14]);

    for (k = 0; k < nln; k = k + 1 )
    {
        for (l = 0; l < nln; l = l + 1 )
        {
            double tmp = 0;
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                tmp = tmp + phi[k+q*nln] * phi[l+q*nln] * w[q];
            }
            LocalMass[k][l] = tmp;
        }
    }

    
    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[4]);

    /* Assembly: loop over the elements */
    int ie;
            
    #pragma omp parallel for shared(invjac,mu,conv_field,si,f,detjac,elements, myRrows, myRcoef,myAcols, myArows, myAcoef, myMcoef) private(gradphi,ie,k,l,q) firstprivate(phi,gradrefphi, w, numRowsElements, nln2, nln, OP, C_t, C_d, LocalMass)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        int d1, d2;
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
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
        }
        
        int iii = 0;
        int ii = 0;
        int a, b;
    
        /* a tes, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            for (b = 0; b < nln; b = b + 1 )
            {
                double aloc = 0;                
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    double diffusion = 0;
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            diffusion = diffusion + C_d[d1][d2] * mu[ie+q*noe] * gradphi[d1][b][q] * gradphi[d2][a][q];
                        }
                    }
                    double transport = 0;
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        transport = transport + C_t[d1] * conv_field[ie+(q+d1*NumQuadPoints)*noe] * gradphi[d1][b][q] * phi[a+q*nln];
                    }
                    
                    double reaction  = si[ie+q*noe] * phi[b+q*nln] * phi[a+q*nln];
                    
                    aloc = aloc + (OP[0] * diffusion + OP[1] * transport + OP[2] * reaction) * w[q];
                }
 
                myArows[ie*nln2+iii] = elements[a+ie*numRowsElements];
                myAcols[ie*nln2+iii] = elements[b+ie*numRowsElements];
                myAcoef[ie*nln2+iii] = aloc*detjac[ie];
                myMcoef[ie*nln2+iii] = LocalMass[a][b]*detjac[ie];
                
                iii = iii + 1;
            }
            
            double floc = 0;
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                floc = floc + ( OP[3] * phi[a+q*nln] * f[ie+q*noe] ) * w[q];
            }
            myRrows[ie*nln+ii] = elements[a+ie*numRowsElements];
            myRcoef[ie*nln+ii] = floc*detjac[ie];
    
            ii = ii + 1;
        }
        
    }
            
}

