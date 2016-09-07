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
#define PRESTRESS_0(ie,q,d1,d2) S_0[ie + q*noe + d1*noe*NumQuadPoints + d2*noe*NumQuadPoints*dim]
#define PRESTRESS_NP1(ie,q,d1,d2) S_np1[ie + q*noe + d1*noe*NumQuadPoints + d2*noe*NumQuadPoints*dim]

#ifndef NHMATERIAL_H_INCLUDED
#define NHMATERIAL_H_INCLUDED

/*************************************************************************/
void NeoHookeanMaterial_forces(mxArray* plhs[], const mxArray* prhs[]);

void NeoHookeanMaterial_jacobian(mxArray* plhs[], const mxArray* prhs[]);

void NeoHookeanMaterial_prestress(mxArray* plhs[], const mxArray* prhs[]);

void NeoHookeanMaterial_jacobianFast(mxArray* plhs[], const mxArray* prhs[]);

void NeoHookeanMaterial_stress(mxArray* plhs[], const mxArray* prhs[]);

#endif

