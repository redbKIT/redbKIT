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
#warning "OpenMP not enabled. Compile with mex CFD_assembler_ExtForces.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=6) {
        mexErrMsgTxt("6 inputs are required.");
    } else if(nlhs>2) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    double* f   = mxGetPr(prhs[0]);
    int noe     = mxGetN(prhs[1]);
    int numRowsElements  = mxGetM(prhs[1]);
    double* elements  = mxGetPr(prhs[1]);

    double* nln_ptr = mxGetPr(prhs[2]);
    int nln     = (int)(nln_ptr[0]);
    
    plhs[0] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
    
    double* myRrows    = mxGetPr(plhs[0]);
    double* myRcoef    = mxGetPr(plhs[1]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[3]);
    
    double* w   = mxGetPr(prhs[3]);
    double* detjac = mxGetPr(prhs[4]);
    double* phi = mxGetPr(prhs[5]);
    
       
    /* Assembly: loop over the elements */
    int ie;
    int a;
    
#pragma omp parallel for shared(detjac,elements,myRrows,myRcoef) private(ie,a,q) firstprivate(phi,w,numRowsElements,nln)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        int ii = 0;
        
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            double floc = 0;
            
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                floc = floc + ( phi[a+q*nln] * f[ie+q*noe] ) * w[q];
            }
            
            myRrows[ie*nln+ii] = elements[a+ie*numRowsElements];
            myRcoef[ie*nln+ii] = floc*detjac[ie];
            ii = ii + 1;
        }
    }
    
}


