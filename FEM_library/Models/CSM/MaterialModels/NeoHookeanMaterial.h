/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */
#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>

#include "Tools.h"

#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]

#ifndef NHMATERIAL_H_INCLUDED
#define NHMATERIAL_H_INCLUDED

/*************************************************************************/
void NeoHookeanMaterial_forces(mxArray* plhs[], const mxArray* prhs[]);

void NeoHookeanMaterial_jacobian(mxArray* plhs[], const mxArray* prhs[]);

void NeoHookeanMaterial_stress(mxArray* plhs[], const mxArray* prhs[]);

#endif

