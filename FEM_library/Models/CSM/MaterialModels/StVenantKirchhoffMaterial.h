/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#include "../../../Core/Tools.h"

#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]

#ifndef SVTMATERIAL_H_INCLUDED
#define SVTMATERIAL_H_INCLUDED

void StVenantKirchhoffMaterial_forces(mxArray* plhs[], const mxArray* prhs[]);

void StVenantKirchhoffMaterial_jacobian(mxArray* plhs[], const mxArray* prhs[]);

void StVenantKirchhoffMaterial_jacobianFast3D(mxArray* plhs[], const mxArray* prhs[]);

void StVenantKirchhoffMaterial_jacobianFast2D(mxArray* plhs[], const mxArray* prhs[]);

void StVenantKirchhoffMaterial_stress(mxArray* plhs[], const mxArray* prhs[]);

#endif