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
#warning "OpenMP not enabled. Compile with mex RBF_evaluate_Fast.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

/*************************************************************************/
double RBF_function(double d, double c, char *RBF_function_name)
{

	double val = -1;

    if (strcmp(RBF_function_name, "gaussian")==0)
	{
		val = exp(-0.5*d*d/(c*c));
        return val;
	}

	if (strcmp(RBF_function_name, "thinplate")==0)
	{
		val = d*d*log(d+1);
        return val;
	}

	if (strcmp(RBF_function_name, "cubic")==0)
	{
		val = (d*d*d);
        return val;
	}
    
    if (strcmp(RBF_function_name, "multiquadric")==0)
	{
		val = sqrt(1+d*d/(c*c));
        return val;
	}

	return val;
     
}
/*************************************************************************/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    char *RBF_function_name = mxArrayToString(prhs[0]);
    
    /* Check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgTxt("5 inputs are required.");
    } else if(nlhs>1) {
        mexErrMsgTxt("Too many output arguments.");
    }

    double* interp_points = mxGetPr(prhs[1]);
    int nI  = mxGetN(prhs[1]);
    
    double* x = mxGetPr(prhs[2]);
    int dimX  = mxGetM(prhs[2]);
    int nPoints  = mxGetN(prhs[2]);
    
    double* tmpPtr = mxGetPr(prhs[3]);
    double constant = tmpPtr[0];
    
    double* coeff = mxGetPr(prhs[4]);
    
    plhs[0] = mxCreateDoubleMatrix(nPoints,1, mxREAL);
    double* I_f    = mxGetPr(plhs[0]);
    
    int i;
    #pragma omp parallel for shared(I_f,x) private(i) firstprivate(coeff,interp_points,nI,dimX,constant,RBF_function_name)
    for (i = 0; i < nPoints; i++)
    {
        int l, k;
        I_f[i] = 0.0;
        for (k = 0; k < nI; k++)
        {
            /*d = distance(x[:,i], interp_points(:,k));*/
            double tmp = 0;
            for (l = 0; l < dimX; l++)
            {
                double tmp2 = (x[l+dimX*i] - interp_points[l+dimX*k]);
                tmp += (tmp2*tmp2);
            }
            double d = sqrt(tmp);
            I_f[i] +=  coeff[k] * RBF_function(d, constant, RBF_function_name);
        }
        
        I_f[i] += coeff[nI];
        
        for (k = 0; k < dimX; k++)
        {
            I_f[i] += coeff[k+nI+1]*x[k+dimX*i];;
        }
    }
    mxFree(RBF_function_name);
}
/*************************************************************************/

