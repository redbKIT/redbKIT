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

#include "LinearElasticMaterial.h"

#ifdef _OPENMP
#include <omp.h>
#else
#warning "OpenMP not enabled. Compile with mex CSM_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

/*************************************************************************/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    
    if(nrhs!=11) {
        mexErrMsgTxt("11 inputs are required.");
    } else if(nlhs>5) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    char *Material_Model = mxArrayToString(prhs[1]);
    
    if (strcmp(Material_Model, "Linear_forces")==0)
    {
            LinearElasticMaterial_forces(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "Linear_jacobian")==0)
    {
            LinearElasticMaterial_jacobian(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "SEMMT_forces")==0)
    {
            SEMMTMaterial_forces(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "SEMMT_jacobian")==0)
    {
            SEMMTMaterial_jacobian(plhs, prhs);
    }
    
   
    /*
    if (strcmp(Material_Model, "NeoHookean_forces")==0)
    {
            NeoHookeanMaterial_forces(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "NeoHookean_jacobian")==0)
    {
            NeoHookeanMaterial_jacobian(plhs, prhs);
    }
    */
    
    if (strcmp(Material_Model, "StVenantKirchhoff_forces")==0)
    {
            StVenantKirchhoffMaterial_forces(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "StVenantKirchhoff_jacobian")==0)
    {
            StVenantKirchhoffMaterial_jacobian(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "NeoHookean2_forces")==0)
    {
            NeoHookean2Material_forces(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "NeoHookean2_jacobian")==0)
    {
            NeoHookean2Material_jacobian(plhs, prhs);
    }
   
    
    /* ===================================== */
    
    if (strcmp(Material_Model, "NeoHookean")==0)
    {
            NeoHookeanMaterial(plhs, prhs);
    }
    
    mxFree(Material_Model);
}
/*************************************************************************/

