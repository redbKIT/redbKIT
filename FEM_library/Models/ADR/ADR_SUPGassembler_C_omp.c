/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch> 
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#include "../../Core/Tools.h"

#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]
#ifdef _OPENMP
    #include <omp.h>
#else
    #warning "OpenMP not enabled. Compile with mex ADR_SUPGassembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=14) {
        mexErrMsgTxt("14 inputs are required.");
    } else if(nlhs>6) {
        mexErrMsgTxt("Too many output arguments.");
    }

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[3]);
    double* nln_ptr = mxGetPr(prhs[4]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[3]);
    int nln2    = nln*nln;
    
    double* tmp_ptr1 = mxGetPr(prhs[2]);
    double dt = tmp_ptr1[0];
    
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
    char *StabType = mxArrayToString(prhs[1]);
    bool flag_t = false;
    if (strcmp(StabType, "SUPGt")==0)
    {
        flag_t = true;
    }
    mxFree(StabType);
    
    /* Local mass matrix (computed only once) with quadrature nodes */
    double LocalMass[nln][nln];
    int q;
    int NumQuadPoints     = mxGetN(prhs[9]);
    
    double* mu   = mxGetPr(prhs[5]);
    double* conv_field   = mxGetPr(prhs[6]);
    double* si   = mxGetPr(prhs[7]);
    double* f    = mxGetPr(prhs[8]);
    double* w   = mxGetPr(prhs[9]);
    double* invjac = mxGetPr(prhs[10]);
    double* detjac = mxGetPr(prhs[11]);
    double* phi = mxGetPr(prhs[12]);
    double* gradrefphi = mxGetPr(prhs[13]);

    int l,k;
    double* elements  = mxGetPr(prhs[3]);

    /* Assembly: loop over the elements */
    int ie;
            
    #pragma omp parallel for shared(invjac,mu,conv_field,si,f,detjac,elements, myRrows, myRcoef,myAcols, myArows, myAcoef, myMcoef) private(ie,k,l,q) firstprivate(phi,gradrefphi, w, numRowsElements, nln2, nln)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        
        double gradphi[dim][nln][NumQuadPoints];

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
        
         /*  compute metric tensors G and g */
        double G[dim][dim];
        double g[dim];
        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            g[d1] = 0;
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                G[d1][d2] = 0.0;
                int d3;
                for (d3 = 0; d3 < dim; d3 = d3 + 1 )
                {
                    G[d1][d2] += INVJAC(ie,d1,d3) * INVJAC(ie,d2,d3);
                }
                g[d1] = g[d1] + INVJAC(ie,d1,d2);
            }
        }
        
        double traceGtG = Mdot(dim, G, G);
        double tauK[NumQuadPoints];
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            double b_hq[dim];
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                b_hq[d1] = conv_field[ie+(q+d1*NumQuadPoints)*noe];
            }
            
            double G_U_hq[dim];
            MatrixVector(dim, dim, G, b_hq, G_U_hq);
             
            tauK[q] = pow( flag_t * 4/(dt*dt) + ScalarProduct(dim, b_hq, G_U_hq) + 9*mu[ie+q*noe]*mu[ie+q*noe]*traceGtG, -0.5);
        }
        
        int iii = 0;
        int ii = 0;
        int a, b;
    
        double bh_gradPHI[nln][NumQuadPoints];
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                bh_gradPHI[k][q] = 0;
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    bh_gradPHI[k][q] += conv_field[ie+(q+d1*NumQuadPoints)*noe] * gradphi[d1][k][q];
                }
            }
        }
        
        /* a test, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            for (b = 0; b < nln; b = b + 1 )
            {
                double aloc = 0;
                double mloc = 0;
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    aloc +=  (bh_gradPHI[b][q] + si[ie+q*noe] * phi[b+q*nln]) * bh_gradPHI[a][q] * tauK[q] * w[q];
                    
                    mloc +=  phi[b+q*nln] * bh_gradPHI[a][q] * tauK[q] * w[q];
                }
 
                myArows[ie*nln2+iii] = elements[a+ie*numRowsElements];
                myAcols[ie*nln2+iii] = elements[b+ie*numRowsElements];
                myAcoef[ie*nln2+iii] = aloc*detjac[ie];
                myMcoef[ie*nln2+iii] = mloc*detjac[ie];
                
                iii = iii + 1;
            }
            
            double floc = 0;
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                floc += ( bh_gradPHI[a][q] * f[ie+q*noe] * tauK[q] ) * w[q];
            }
            myRrows[ie*nln+ii] = elements[a+ie*numRowsElements];
            myRcoef[ie*nln+ii] = floc*detjac[ie];
    
            ii = ii + 1;
        }
        
    }
            
}

