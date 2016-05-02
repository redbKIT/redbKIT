/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch> 
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#ifdef _OPENMP
    #include <omp.h>
#else
    #warning "OpenMP not enabled. Compile with mex ADR_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=6) {
        mexErrMsgTxt("6 inputs are required.");
    } else if(nlhs>3) {
        mexErrMsgTxt("Too many output arguments.");
    }

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[1]);
    double* nln_ptr = mxGetPr(prhs[2]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[1]);
    int nln2    = nln*nln;
    
    /**/
    plhs[0] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
    
    double* myMrows    = mxGetPr(plhs[0]);
    double* myMcols    = mxGetPr(plhs[1]);
    double* myMcoef    = mxGetPr(plhs[2]);
        
    /* Local mass matrix (computed only once) with quadrature nodes */
    double LocalMass[nln][nln];
    int q;
    int NumQuadPoints     = mxGetN(prhs[3]);
    
    double* w   = mxGetPr(prhs[3]);
    double* detjac = mxGetPr(prhs[4]);
    double* phi = mxGetPr(prhs[5]);

    int k, l;
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

    double* elements  = mxGetPr(prhs[1]);

    /* Assembly: loop over the elements */
    int ie;
            
    #pragma omp parallel for shared(detjac,elements,myMcols, myMrows, myMcoef) private(ie) firstprivate(phi, numRowsElements, nln2, LocalMass)
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
       
        int iii = 0;
        int a, b;
    
        /* a tes, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            for (b = 0; b < nln; b = b + 1 )
            {
                myMrows[ie*nln2+iii] = elements[a+ie*numRowsElements];
                myMcols[ie*nln2+iii] = elements[b+ie*numRowsElements];
                myMcoef[ie*nln2+iii] = LocalMass[a][b]*detjac[ie];
                
                iii = iii + 1;
            }
        }
        
    }
            
}

