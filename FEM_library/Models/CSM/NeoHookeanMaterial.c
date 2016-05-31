/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "NeoHookeanMaterial.h"

/*************************************************************************/
void NeoHookeanMaterial(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    double* myRrows    = mxGetPr(plhs[3]);
    double* myRcoef    = mxGetPr(plhs[4]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);
    
    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);
    
    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[4]);
    
    double GradV[dim][dim];
    double GradU[dim][dim];
    double GradUh[dim][dim][NumQuadPoints];
    
    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }
    
    double F[NumQuadPoints][dim][dim];
    double invFT[NumQuadPoints][dim][dim];
    double detF[NumQuadPoints];
    double logdetF[NumQuadPoints];
    double dP[dim][dim];
    double P_Uh[dim][dim];
    
    double dF[dim][dim];
    double dE[dim][dim];
    
    double* material_param = mxGetPr(prhs[2]);
    double Young = material_param[0];
    double Poisson = material_param[1];
    double mu = Young / (2 + 2 * Poisson);
    double lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myRrows,myRcoef,myAcols,myArows,myAcoef,U_h) private(gradphi,F,invFT,logdetF,detF,dP,P_Uh,dF,GradV,GradU,GradUh,ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {             
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nln; k = k + 1 )
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
            
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[d1][d2][q] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[d1][d2][q] = GradUh[d1][d2][q] + U_h[e_k] * gradphi[d2][k][q];
                    }
                    F[q][d1][d2]  = Id[d1][d2] + GradUh[d1][d2][q];
                }
            }
            detF[q] = MatrixDeterminant(dim, F[q]);
            MatrixInvT(dim, F[q], invFT[q] );
            logdetF[q] = log( detF[q] );
        }

        int iii = 0;
        int ii = 0;
        int a, b, i_c, j_c;
        
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
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
                
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
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
                                GradV[i_c][d2] = gradphi[d2][a][q];
                                GradU[j_c][d2] = gradphi[d2][b][q];
                            }
                            
                            
                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    dF[d1][d2] = GradU[d1][d2];
                                }
                            }
                                                        
                            double P2[dim][dim];
                            double P3[dim][dim];
                            double P_tmp[dim][dim];
                            double P_tmp2[dim][dim];
                            
                            MatrixProductAlphaT2(dim, mu-lambda*logdetF[q], invFT[q], dF, P_tmp);
                            MatrixProductAlpha(dim, 1.0, P_tmp, invFT[q], dP);
                            
                            MatrixProductAlphaT1( dim, 1.0, invFT[q], dF, P_tmp2);
                            double coef = lambda * Trace(dim, P_tmp2);
                            
                            MatrixProductAlpha(dim, coef, Id, invFT[q], P2 );
                            
                            MatrixSum(dim, dP, P2);
                            
                            MatrixProductAlpha(dim, mu, Id, dF, P3 );
                            
                            MatrixSum(dim, dP, P3);
                            
                            aloc  = aloc + Mdot( dim, GradV, dP) * w[q];
                        }
                        myArows[ie*nln2*dim*dim+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*dim*dim+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*dim*dim+iii] = aloc*detjac[ie];
                        
                        iii = iii + 1;
                    }
                }
                
                double rloc = 0;
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[i_c][d2] = gradphi[d2][a][q];
                    }
                    
                    double P1[dim][dim];
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            P_Uh[d1][d2] =  ( mu * ( F[q][d1][d2] - invFT[q][d1][d2] ) + lambda * logdetF[q] * invFT[q][d1][d2] );
                        }
                    }
                    rloc  = rloc + Mdot( dim, GradV, P_Uh) * w[q];
                }
                                            
                myRrows[ie*nln*dim+ii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                myRcoef[ie*nln*dim+ii] = rloc*detjac[ie];
                ii = ii + 1;
            }
        }
    }
        
}
/*************************************************************************/

