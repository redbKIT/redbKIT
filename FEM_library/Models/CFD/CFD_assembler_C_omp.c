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
#define GRADREFPHIV(i,j,k) gradrefphiV[i+(j+k*NumQuadPoints)*nlnV]
#ifdef _OPENMP
#include <omp.h>
#else
#warning "OpenMP not enabled. Compile with mex CFD_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
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
void AssembleStokes(mxArray* plhs[], const mxArray* prhs[])
{
    double* dim_ptr = mxGetPr(prhs[2]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[3]);
    double* nln_ptrV = mxGetPr(prhs[4]);
    int nlnV     = (int)(nln_ptrV[0]);
    double* nln_ptrP = mxGetPr(prhs[5]);
    int nlnP     = (int)(nln_ptrP[0]);
    int numRowsElements  = mxGetM(prhs[3]);
    
    int nln  = nlnV + nlnP;
    int nln2 = nln*nln;
    int local_matrix_size = nlnV*nlnV*dim*dim + 2*nlnV*nlnP*dim;
    int global_lenght = noe * local_matrix_size;
    
    plhs[0] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[7]);
    
    double* NumNodes_ptr = mxGetPr(prhs[6]);
    int NumScalarDofsV     = (int)(NumNodes_ptr[0] / dim);
    
    double* w   = mxGetPr(prhs[7]);
    double* invjac = mxGetPr(prhs[8]);
    double* detjac = mxGetPr(prhs[9]);
    double* phiV = mxGetPr(prhs[10]);
    double* gradrefphiV = mxGetPr(prhs[11]);
    double* phiP = mxGetPr(prhs[12]);
    
    double gradphiV[dim][nlnV][NumQuadPoints];
    double* elements  = mxGetPr(prhs[3]);
    
    double GradV[dim][dim];
    double GradU[dim][dim];
    double F[dim][dim];
    
    double* material_param = mxGetPr(prhs[1]);
    double viscosity = material_param[0];
    
    /* Assembly: loop over the elements */
    int ie, d1, d2;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef) private(gradphiV,F,GradV,GradU,ie,k,l,q,d1,d2) firstprivate(phiV,phiP,gradrefphiV,w,numRowsElements,local_matrix_size,nlnV,nlnP,NumScalarDofsV,viscosity)
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        for (k = 0; k < nlnV; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphiV[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphiV[d1][k][q] = gradphiV[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHIV(k,q,d2);
                    }
                }
            }
        }
        
        int iii = 0;
        int ii = 0;
        int a, b, i_c, j_c;
        
        /* loop over velocity test functions --> a */
        for (a = 0; a < nlnV; a = a + 1 )
        {
            /* loop over test velocity components --> i_c */
            for (i_c = 0; i_c < dim; i_c = i_c + 1 )
            {
                /* set gradV to zero*/
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[d1][d2] = 0;
                    }
                }
                
                /* loop over velocity trial functions --> b */
                for (b = 0; b < nlnV; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < dim; j_c = j_c + 1 )
                    {
                        /* set gradU to zero*/
                        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                        {
                            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                            {
                                GradU[d1][d2] = 0;
                            }
                        }
                        
                        double aloc = 0;
                        for (q = 0; q < NumQuadPoints; q = q + 1 )
                        {
                            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                            {
                                GradV[i_c][d2] = gradphiV[d2][a][q];
                                GradU[j_c][d2] = gradphiV[d2][b][q];
                            }
                            
                            
                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    F[d1][d2] = viscosity *  ( GradU[d1][d2] + GradU[d2][d1] );
                                }
                            }
                            
                            aloc  = aloc + Mdot( dim, GradV, F) * w[q];
                        }
                        myArows[ie*local_matrix_size+iii] = elements[a+ie*numRowsElements] + i_c * NumScalarDofsV;
                        myAcols[ie*local_matrix_size+iii] = elements[b+ie*numRowsElements] + j_c * NumScalarDofsV;
                        myAcoef[ie*local_matrix_size+iii] = aloc*detjac[ie];
            
                        iii = iii + 1;
                    }
                }
                
                /* set gradV to zero*/
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[d1][d2] = 0;
                    }
                }
                
                
                /* loop over pressure trial functions --> b */
                for (b = 0; b < nlnP; b = b + 1 )
                {
                    double aloc = 0;
                    for (q = 0; q < NumQuadPoints; q = q + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            GradV[i_c][d2] = gradphiV[d2][a][q];
                        }
                        
                        aloc  = aloc + ( - phiP[b+q*nlnP] * Trace(dim, GradV) ) * w[q];
                    }
                    
                    
                    myArows[ie*local_matrix_size+iii] = elements[a+ie*numRowsElements] + i_c * NumScalarDofsV;
                    myAcols[ie*local_matrix_size+iii] = elements[b+ie*numRowsElements] + dim * NumScalarDofsV;
                    myAcoef[ie*local_matrix_size+iii] = aloc*detjac[ie];
                    
                    iii = iii + 1;
                    
                    myAcols[ie*local_matrix_size+iii] = elements[a+ie*numRowsElements] + i_c * NumScalarDofsV;
                    myArows[ie*local_matrix_size+iii] = elements[b+ie*numRowsElements] + dim * NumScalarDofsV;
                    myAcoef[ie*local_matrix_size+iii] = -aloc*detjac[ie];
                    
                    iii = iii + 1;
                    
                    
                }
            }
        }
    }
    
}

/*************************************************************************/
void AssembleConvective_Oseen(mxArray* plhs[], const mxArray* prhs[])
{
    double* dim_ptr = mxGetPr(prhs[2]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[3]);
    double* nln_ptrV = mxGetPr(prhs[4]);
    int nlnV     = (int)(nln_ptrV[0]);
    int numRowsElements  = mxGetM(prhs[3]);
    
    int local_matrix_size = nlnV*nlnV*dim;
    int global_lenght = noe * local_matrix_size;
    
    plhs[0] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    
    double* NumNodes_ptr = mxGetPr(prhs[5]);
    int NumScalarDofsV     = (int)(NumNodes_ptr[0] / dim);
    
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phiV = mxGetPr(prhs[9]);
    double* gradrefphiV = mxGetPr(prhs[10]);
    double* U_h   = mxGetPr(prhs[11]);
    
    double gradphiV[NumQuadPoints][dim][nlnV];
    double* elements  = mxGetPr(prhs[3]);
    
    double GradV[dim][dim];
    double GradU[dim][dim];
    double U_hq[NumQuadPoints][dim];
    
    double* material_param = mxGetPr(prhs[1]);
    double density = material_param[0];
    
    /* Assembly: loop over the elements */
    int ie, d1, d2;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(gradphiV,GradV,GradU,U_hq,ie,k,l,q,d1,d2) firstprivate(phiV,gradrefphiV,w,numRowsElements,local_matrix_size,nlnV,NumScalarDofsV,density)
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (k = 0; k < nlnV; k = k + 1 )
                {
                    gradphiV[q][d1][k] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphiV[q][d1][k] = gradphiV[q][d1][k] + INVJAC(ie,d1,d2)*GRADREFPHIV(k,q,d2);
                    }
                }
            }
        }

        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                U_hq[q][d1] = 0;
                for (k = 0; k < nlnV; k = k + 1 )
                {
                    int e_k;
                    e_k = (int)(elements[ie*numRowsElements + k] + d1*NumScalarDofsV - 1);
                    U_hq[q][d1] = U_hq[q][d1] + U_h[e_k] * phiV[k+q*nlnV];
                }
            }
        }
        
        int iii = 0;
        int a, b, i_c, j_c;
        
        /* loop over velocity test functions --> a */
        for (a = 0; a < nlnV; a = a + 1 )
        {
            /* loop over velocity trial functions --> b */
            for (b = 0; b < nlnV; b = b + 1 )
            {
                double aloc = 0;
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        aloc  = aloc + U_hq[q][d1] * gradphiV[q][d1][b] * phiV[a+q*nlnV] * w[q];
                    }
                }
                
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    myArows[ie*local_matrix_size+iii] = elements[a+ie*numRowsElements] + d1 * NumScalarDofsV;
                    myAcols[ie*local_matrix_size+iii] = elements[b+ie*numRowsElements] + d1 * NumScalarDofsV;
                    myAcoef[ie*local_matrix_size+iii] = density*aloc*detjac[ie];
                    
                    iii = iii + 1;
                }
            }
        }
    }
    
}
/*************************************************************************/
void AssembleConvective(mxArray* plhs[], const mxArray* prhs[])
{
    double* dim_ptr = mxGetPr(prhs[2]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[3]);
    double* nln_ptrV = mxGetPr(prhs[4]);
    int nlnV     = (int)(nln_ptrV[0]);
    int numRowsElements  = mxGetM(prhs[3]);
    
    int local_matrix_size1 = nlnV*nlnV*dim;
    int local_matrix_size2 = nlnV*nlnV*dim*dim;

    int global_lenght1 = noe * local_matrix_size1;
    int global_lenght2 = noe * local_matrix_size2;

    plhs[0] = mxCreateDoubleMatrix(global_lenght1,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(global_lenght1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(global_lenght1,1, mxREAL);
    
    plhs[3] = mxCreateDoubleMatrix(global_lenght2,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(global_lenght2,1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(global_lenght2,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    double* myBrows    = mxGetPr(plhs[3]);
    double* myBcols    = mxGetPr(plhs[4]);
    double* myBcoef    = mxGetPr(plhs[5]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    
    double* NumNodes_ptr = mxGetPr(prhs[5]);
    int NumScalarDofsV     = (int)(NumNodes_ptr[0] / dim);
    
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phiV = mxGetPr(prhs[9]);
    double* gradrefphiV = mxGetPr(prhs[10]);
    double* U_h   = mxGetPr(prhs[11]);
    
    double gradphiV[dim][nlnV][NumQuadPoints];
    double* elements  = mxGetPr(prhs[3]);
    
    double GradV[dim][dim];
    double GradU[dim][dim];
    double U_hq[dim][NumQuadPoints];
    double GradUh[dim][dim][NumQuadPoints];
    
    double* material_param = mxGetPr(prhs[1]);
    double density = material_param[0];
    
    /* Assembly: loop over the elements */
    int ie, d1, d2;
        
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(gradphiV,GradV,GradU,GradUh,U_hq,ie,k,l,q,d1,d2) firstprivate(phiV,gradrefphiV,w,numRowsElements,local_matrix_size1,local_matrix_size2,nlnV,NumScalarDofsV,density)
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nlnV; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphiV[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphiV[d1][k][q] = gradphiV[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHIV(k,q,d2);
                    }
                }
            }
            /* Compute U_h and Grad(U_h) on the quadrature nodes of the current element*/
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                U_hq[d1][q] = 0;
                for (k = 0; k < nlnV; k = k + 1 )
                {
                    int e_k;
                    e_k = (int)(elements[ie*numRowsElements + k] + d1*NumScalarDofsV - 1);
                    U_hq[d1][q] = U_hq[d1][q] + U_h[e_k] * phiV[k+q*nlnV];
                }
                
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[d1][d2][q] = 0;
                    for (k = 0; k < nlnV; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumScalarDofsV - 1);
                        GradUh[d1][d2][q] = GradUh[d1][d2][q] + U_h[e_k] * gradphiV[d2][k][q];
                    }
                }
            }
        }
        
        int iii = 0, iii2 = 0;
        int a, b, i_c, j_c;
        
        /* loop over velocity test functions --> a */
        for (a = 0; a < nlnV; a = a + 1 )
        {
            /* Assemble C1 = U_h * grad(u) * v */
            /* loop over velocity trial functions --> b */
            for (b = 0; b < nlnV; b = b + 1 )
            {
                double aloc = 0;
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        aloc  = aloc + U_hq[d1][q] * gradphiV[d1][b][q] * phiV[a+q*nlnV] * w[q];
                    }
                }
                
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    myArows[ie*local_matrix_size1+iii] = elements[a+ie*numRowsElements] + d1 * NumScalarDofsV;
                    myAcols[ie*local_matrix_size1+iii] = elements[b+ie*numRowsElements] + d1 * NumScalarDofsV;
                    myAcoef[ie*local_matrix_size1+iii] = density*aloc*detjac[ie];
                    
                    iii = iii + 1;
                }
            }
            
            /* Assemble C2 = u * grad(U_h) * v */
            /* loop over test velocity components --> i_c */
            for (i_c = 0; i_c < dim; i_c = i_c + 1 )
            {
                /* loop over velocity trial functions --> b */
                for (b = 0; b < nlnV; b = b + 1 )
                {
                    /* loop over test velocity components --> j_c */
                    for (j_c = 0; j_c < dim; j_c = j_c + 1 )
                    {
                        double aloc = 0;
                        for (q = 0; q < NumQuadPoints; q = q + 1 )
                        {
                            aloc  = aloc + phiV[b+q*nlnV] * GradUh[i_c][j_c][q] * phiV[a+q*nlnV] * w[q];
                        }
                        myBrows[ie*local_matrix_size2+iii2] = elements[a+ie*numRowsElements] + i_c * NumScalarDofsV;
                        myBcols[ie*local_matrix_size2+iii2] = elements[b+ie*numRowsElements] + j_c * NumScalarDofsV;
                        myBcoef[ie*local_matrix_size2+iii2] = density*aloc*detjac[ie];
                        
                        iii2 = iii2 + 1;
                        
                    }
                }
            }

                
        }
    }
}
/*************************************************************************/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    char *Assembly_name = mxArrayToString(prhs[0]);
    
    if (strcmp(Assembly_name, "Stokes")==0)
    {
        /* Check for proper number of arguments */
        if(nrhs!=13) {
            mexErrMsgTxt("13 inputs are required.");
        } else if(nlhs>3) {
            mexErrMsgTxt("Too many output arguments.");
        }

        AssembleStokes(plhs, prhs);
    }
    
    
    if (strcmp(Assembly_name, "convective_Oseen")==0)
    {
        /* Check for proper number of arguments */
        if(nrhs!=12) {
            mexErrMsgTxt("12 inputs are required.");
        } else if(nlhs>3) {
            mexErrMsgTxt("Too many output arguments.");
        }
        
        AssembleConvective_Oseen(plhs, prhs);
    }
    if (strcmp(Assembly_name, "convective")==0)
    {
        /* Check for proper number of arguments */
        if(nrhs!=12) {
            mexErrMsgTxt("12 inputs are required.");
        } else if(nlhs>6) {
            mexErrMsgTxt("Too many output arguments.");
        }
        
        AssembleConvective(plhs, prhs);
    }
    
    mxFree(Assembly_name);
}
/*************************************************************************/

