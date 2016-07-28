/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "RaghavanVorpMaterial.h"

/*************************************************************************/
void RaghavanVorpMaterial_forces(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    
    double* myRrows    = mxGetPr(plhs[0]);
    double* myRcoef    = mxGetPr(plhs[1]);
    
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
    
    double* elements  = mxGetPr(prhs[4]);
    
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
    
    double* material_param = mxGetPr(prhs[2]);
    double alpha = material_param[0];
    double beta = material_param[1];
    double bulk = material_param[2];
    
    /*
    double mu = Young / (2.0 + 2.0 * Poisson);
    double lambda =  Young * Poisson /( (1.0 + Poisson) * (1.0-2.0*Poisson) );
    double bulk = ( 2.0 / 3.0 ) * mu + lambda;
    */
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myRrows,myRcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,NumQuadPoints,numRowsElements,nln2,nln,NumNodes,Id,alpha,beta,bulk,noe,dim)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {             
        double I_C[NumQuadPoints];
        double detF[NumQuadPoints];
        double logdetF[NumQuadPoints];
        double pow23detF[NumQuadPoints];
        double pow2detF[NumQuadPoints];
        
        double F[NumQuadPoints][dim][dim];
        double invFT[NumQuadPoints][dim][dim];
        double C[NumQuadPoints][dim][dim];
        
        double dP[dim][dim];
        double P_Uh[dim][dim];
        
        double GradV[dim][dim];
        double GradUh[NumQuadPoints][dim][dim];
        
        double gradphi[dim][nln][NumQuadPoints];

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
                    GradUh[q][d1][d2] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[q][d1][d2] = GradUh[q][d1][d2] + U_h[e_k] * gradphi[d2][k][q];
                    }
                    F[q][d1][d2]  = Id[d1][d2] + GradUh[q][d1][d2];
                }
            }
            detF[q] = MatrixDeterminant(dim, F[q]);
            MatrixInvT(dim, F[q], invFT[q] );
            MatrixProductAlphaT1(dim, 1.0, F[q], F[q], C[q] );
            logdetF[q] = log( detF[q] );
            pow23detF[q] = pow(detF[q], -2.0 / 3.0);
            pow2detF[q] = pow(detF[q], 2.0);
            I_C[q] = Trace(dim, C[q]);
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
                            P_Uh[d1][d2] =  2.0 * ( alpha + 2.0 * beta * ( pow23detF[q] * I_C[q] - 3.0 ) )   
                                              * pow23detF[q] * ( F[q][d1][d2] - 1.0 / 3.0 * I_C[q]  * invFT[q][d1][d2] ) 
                                            + 1.0 / 2.0 * bulk * ( pow2detF[q] - detF[q] + logdetF[q] ) * invFT[q][d1][d2];                      
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


/*************************************************************************/
void RaghavanVorpMaterial_jacobian(mxArray* plhs[], const mxArray* prhs[])
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
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
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
    
    double* elements  = mxGetPr(prhs[4]);
    
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
    
    double* material_param = mxGetPr(prhs[2]);
    double alpha = material_param[0];
    double beta = material_param[1];
    double bulk = material_param[2];
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,alpha,beta,bulk)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {             
        double I_C[NumQuadPoints];
        double detF[NumQuadPoints];
        double logdetF[NumQuadPoints];
        double pow23detF[NumQuadPoints];
        double pow2detF[NumQuadPoints];
        
        double F[NumQuadPoints][dim][dim];
        double invFT[NumQuadPoints][dim][dim];
        double C[NumQuadPoints][dim][dim];
        
        double dP[dim][dim];
        double P_Uh[dim][dim];
        
        double GradV[dim][dim];
        double GradU[dim][dim];
        double GradUh[NumQuadPoints][dim][dim];
        
        double gradphi[dim][nln][NumQuadPoints];

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
                    GradUh[q][d1][d2] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[q][d1][d2] = GradUh[q][d1][d2] + U_h[e_k] * gradphi[d2][k][q];
                    }
                    F[q][d1][d2]  = Id[d1][d2] + GradUh[q][d1][d2];
                }
            }
            detF[q] = MatrixDeterminant(dim, F[q]);
            MatrixInvT(dim, F[q], invFT[q] );
            MatrixProductAlphaT1(dim, 1.0, F[q], F[q], C[q] );
            logdetF[q] = log( detF[q] );
            pow23detF[q] = pow(detF[q], -2.0 / 3.0);
            pow2detF[q] = pow(detF[q], 2.0);
            I_C[q] = Trace(dim, C[q]);
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
                            
                            /* volumetric part */
                            double dP_vol[dim][dim];
                            double dP_vol1[dim][dim];
                            double dP_vol2_tmp[dim][dim];
                            double dP_vol2[dim][dim];
                            
                            MatrixScalar(dim, 0.5*bulk * (2.0*pow2detF[q] -detF[q] + 1.0)*Mdot(dim, invFT[q], GradU),
                                    invFT[q], dP_vol);
                            
                            MatrixProductAlphaT2(dim, 0.5*bulk * ( - pow2detF[q] + detF[q] - logdetF[q]), invFT[q], GradU, dP_vol2_tmp);
                            MatrixProductAlpha(dim, 1.0, dP_vol2_tmp, invFT[q], dP_vol2);
                            
                            MatrixSum(dim, dP_vol, dP_vol2);
                            
                            /* isochoric part */
                            double dP_iso[dim][dim];
                            double dP_iso1[dim][dim];
                            double dP_iso24[dim][dim];
                            double dP_iso3[dim][dim];
                            double dP_iso5[dim][dim];
                            double dP_iso5_tmp[dim][dim];
                            double dP_iso5_tmp2[dim][dim];
                            
                            double mu_q = 2.0 * ( alpha + 2.0 * beta * ( pow23detF[q] * I_C[q] - 3.0 ) );
                            
                            MatrixScalar(dim, -2.0 / 3.0 * mu_q * pow23detF[q] * Mdot(dim, invFT[q], GradU),
                                    F[q], dP_iso1);
                            
                            MatrixScalar(dim, mu_q * pow23detF[q] * 
                                                ( 2.0 / 9.0 * I_C[q]  * Mdot(dim, invFT[q], GradU) 
                                                 -2.0 / 3.0 * Mdot(dim, F[q], GradU) ),
                                    invFT[q], dP_iso24);
                            
                            MatrixScalar(dim, mu_q * pow23detF[q], GradU, dP_iso3);
                            
                            
                            MatrixProductAlphaT2(dim, 1.0, invFT[q], GradU, dP_iso5_tmp);
                            MatrixProductAlpha(dim, 1.0, dP_iso5_tmp, invFT[q], dP_iso5_tmp2);
                            MatrixScalar(dim, 1.0 / 3.0 * mu_q * pow23detF[q] * I_C[q] , dP_iso5_tmp2, dP_iso5);
                                                        
                            /* multiplicative factor: */
                            double dP_iso6[dim][dim];
                            
                            double dP_6_tmp[dim][dim];
                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    dP_6_tmp[d1][d2] =   F[q][d1][d2] - 1.0 / 3.0 * I_C[q]  * invFT[q][d1][d2] ;
                                }
                            }

                            double scalar = 2.0 *  pow23detF[q] * 2.0 * beta * 2.0 * pow23detF[q] * Mdot(dim, dP_6_tmp,  GradU);
                            MatrixScalar(dim, scalar , dP_6_tmp, dP_iso6);

                            /* Sum all contributes */
                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    dP[d1][d2] = dP_vol[d1][d2] 
                                                + dP_iso1[d1][d2] 
                                                + dP_iso24[d1][d2]
                                                + dP_iso3[d1][d2] 
                                                + dP_iso5[d1][d2]
                                                + dP_iso6[d1][d2];
                                }
                            }
                            aloc  = aloc + Mdot( dim, GradV, dP) * w[q];
                        }
                        myArows[ie*nln2*dim*dim+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*dim*dim+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*dim*dim+iii] = aloc*detjac[ie];
                        
                        iii = iii + 1;
                    }
                }
                
            }
        }
    }
        
}

/*************************************************************************/
void RaghavanVorpMaterial_jacobianFast(mxArray* plhs[], const mxArray* prhs[])
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
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
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
    
    double* elements  = mxGetPr(prhs[4]);
    
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
    
    double* material_param = mxGetPr(prhs[2]);
    double alpha = material_param[0];
    double beta = material_param[1];
    double bulk = material_param[2];
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,alpha,beta,bulk)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {             
        double I_C[NumQuadPoints];
        double detF[NumQuadPoints];
        double logdetF[NumQuadPoints];
        double pow23detF[NumQuadPoints];
        double pow2detF[NumQuadPoints];
        
        double F[NumQuadPoints][dim][dim];
        double invFT[NumQuadPoints][dim][dim];
        double C[NumQuadPoints][dim][dim];
        
        double GradUh[NumQuadPoints][dim][dim];
        double gradphi[NumQuadPoints][dim][nln];

        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nln; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[q][d1][k] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[q][d1][k] = gradphi[q][d1][k] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
            
            for (d1 = 0; d1 < 3; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < 3; d2 = d2 + 1 )
                {
                    GradUh[q][d1][d2] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[q][d1][d2] = GradUh[q][d1][d2] + U_h[e_k] * gradphi[q][d2][k];
                    }
                    F[q][d1][d2]  = Id[d1][d2] + GradUh[q][d1][d2];
                }
            }
            detF[q] = MatrixDeterminant3(dim, F[q]);
            MatrixInvT3(dim, F[q], invFT[q] );
            MatrixProductAlphaT1(dim, 1.0, F[q], F[q], C[q] );
            logdetF[q] = log( detF[q] );
            pow23detF[q] = pow(detF[q], -2.0 / 3.0);
            pow2detF[q] = pow(detF[q], 2.0);
            I_C[q] = Trace(dim, C[q]);
        }

        int iii = 0;
        int a, b, i_c, j_c;
        
        double aloc[nln][dim][nln][dim];
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 3; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 3; j_c = j_c + 1 )
                    {
                        aloc[a][i_c][b][j_c]  = 0.0;
                    }
                }
            }
        }
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            double mu_q = 2.0 * ( alpha + 2.0 * beta * ( pow23detF[q] * I_C[q] - 3.0 ) );
            
            double vol_factor1 = 0.5*bulk * (2.0*pow2detF[q] -detF[q] + 1.0);
            double vol_factor2 = 0.5*bulk * ( - pow2detF[q] + detF[q] - logdetF[q]);
                    
            double P_F[dim][dim];
            
            for (d1 = 0; d1 < 3; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < 3; d2 = d2 + 1 )
                {
                    P_F[d1][d2] =   F[q][d1][d2] - 1.0 / 3.0 * I_C[q]  * invFT[q][d1][d2] ;
                }
            }

            /* loop over test functions --> a */
            for (a = 0; a < nln; a = a + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {

                    aloc[a][0][b][0] += ( gradphi[q][0][a]*(invFT[q][0][0]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][0]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][0]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][0][b] - (2*P_F[0][0]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][0]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][0][1]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][1]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][1]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][1][b] - (2*P_F[0][1]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][1]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][0][2]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][2]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][2]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][2][b] - (2*P_F[0][2]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][2]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) ) * w[q];
                    
                    aloc[a][0][b][1] += ( gradphi[q][0][a]*(invFT[q][1][0]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][0]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][0]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) - (2*P_F[0][0]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][0]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][1][1]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][1]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][1]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) - (2*P_F[0][1]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][1]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][1][2]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][2]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][2]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) - (2*P_F[0][2]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][2]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) ) * w[q];
                    
                    aloc[a][0][b][2] += ( gradphi[q][0][a]*(invFT[q][2][0]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][0]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][0]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) - (2*P_F[0][0]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][0]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][2][1]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][1]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][1]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) - (2*P_F[0][1]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][1]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][2][2]*(invFT[q][0][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][0][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][0][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][2]*((I_C[q]*invFT[q][0][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][0][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][0][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][0][2]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) - (2*P_F[0][2]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][0][2]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[0][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) ) * w[q];
                    
                    aloc[a][1][b][0] += ( gradphi[q][0][a]*(invFT[q][0][0]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][0]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][0]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) - (2*P_F[1][0]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][0]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][0][1]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][1]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][1]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) - (2*P_F[1][1]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][1]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][0][2]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][2]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][2]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) - (2*P_F[1][2]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][2]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) ) * w[q];
                    
                    aloc[a][1][b][1] += ( gradphi[q][0][a]*(invFT[q][1][0]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][0]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][0]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][0][b] - (2*P_F[1][0]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][0]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][1][1]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][1]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][1]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][1][b] - (2*P_F[1][1]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][1]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][1][2]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][2]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][2]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][2][b] - (2*P_F[1][2]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][2]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) ) * w[q];
                    
                    aloc[a][1][b][2] += ( gradphi[q][0][a]*(invFT[q][2][0]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][0]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][0]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) - (2*P_F[1][0]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][0]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][2][1]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][1]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][1]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) - (2*P_F[1][1]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][1]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][2][2]*(invFT[q][1][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][1][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][1][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][2]*((I_C[q]*invFT[q][1][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][1][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][1][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][1][2]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) - (2*P_F[1][2]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][1][2]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[1][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) ) * w[q];
                    
                    aloc[a][2][b][0] += ( gradphi[q][0][a]*(invFT[q][0][0]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][0]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][0]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) - (2*P_F[2][0]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][0]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][0][1]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][1]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][1]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) - (2*P_F[2][1]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][1]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][0][2]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][0][2]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][2]*vol_factor1*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]) - (2*P_F[2][2]*mu_q*pow23detF[q]*(invFT[q][0][0]*gradphi[q][0][b] + invFT[q][0][1]*gradphi[q][1][b] + invFT[q][0][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][2]*mu_q*pow23detF[q]*(F[q][0][0]*gradphi[q][0][b] + F[q][0][1]*gradphi[q][1][b] + F[q][0][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[0][0]*gradphi[q][0][b] + P_F[0][1]*gradphi[q][1][b] + P_F[0][2]*gradphi[q][2][b])) ) * w[q];
                    
                    aloc[a][2][b][1] += ( gradphi[q][0][a]*(invFT[q][1][0]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][0]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][0]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) - (2*P_F[2][0]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][0]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][1][1]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][1]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][1]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) - (2*P_F[2][1]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][1]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][1][2]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][1][2]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][2]*vol_factor1*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]) - (2*P_F[2][2]*mu_q*pow23detF[q]*(invFT[q][1][0]*gradphi[q][0][b] + invFT[q][1][1]*gradphi[q][1][b] + invFT[q][1][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][2]*mu_q*pow23detF[q]*(F[q][1][0]*gradphi[q][0][b] + F[q][1][1]*gradphi[q][1][b] + F[q][1][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[1][0]*gradphi[q][0][b] + P_F[1][1]*gradphi[q][1][b] + P_F[1][2]*gradphi[q][2][b])) ) * w[q];
                    
                    aloc[a][2][b][2] += ( gradphi[q][0][a]*(invFT[q][2][0]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][0]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][0]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][0][b] - (2*P_F[2][0]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][0]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][0]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) + gradphi[q][1][a]*(invFT[q][2][1]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][1]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][1]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][1][b] - (2*P_F[2][1]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][1]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][1]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) + gradphi[q][2][a]*(invFT[q][2][2]*(invFT[q][2][0]*gradphi[q][0][b]*vol_factor2 + invFT[q][2][1]*gradphi[q][1][b]*vol_factor2 + invFT[q][2][2]*gradphi[q][2][b]*vol_factor2) + invFT[q][2][2]*((I_C[q]*invFT[q][2][0]*mu_q*pow23detF[q]*gradphi[q][0][b])/3.0 + (I_C[q]*invFT[q][2][1]*mu_q*pow23detF[q]*gradphi[q][1][b])/3.0 + (I_C[q]*invFT[q][2][2]*mu_q*pow23detF[q]*gradphi[q][2][b])/3.0) + invFT[q][2][2]*vol_factor1*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]) + mu_q*pow23detF[q]*gradphi[q][2][b] - (2*P_F[2][2]*mu_q*pow23detF[q]*(invFT[q][2][0]*gradphi[q][0][b] + invFT[q][2][1]*gradphi[q][1][b] + invFT[q][2][2]*gradphi[q][2][b]))/3.0 - (2*invFT[q][2][2]*mu_q*pow23detF[q]*(F[q][2][0]*gradphi[q][0][b] + F[q][2][1]*gradphi[q][1][b] + F[q][2][2]*gradphi[q][2][b]))/3.0 + 8*P_F[2][2]*beta*pow23detF[q]*pow23detF[q]*(P_F[2][0]*gradphi[q][0][b] + P_F[2][1]*gradphi[q][1][b] + P_F[2][2]*gradphi[q][2][b])) ) * w[q];

                }
            }
        }

        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 3; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 3; j_c = j_c + 1 )
                    {
                        myArows[ie*nln2*9+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*9+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*9+iii] = aloc[a][i_c][b][j_c]*detjac[ie];
                        iii = iii + 1;
                    }
                }
            }
        }
    }
        
}
/*************************************************************************/

void RaghavanVorpMaterial_stress(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(noe,dim*dim, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(noe,dim*dim, mxREAL);
    
    double* P        = mxGetPr(plhs[0]);
    double* Sigma    = mxGetPr(plhs[1]);

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
    
    double* material_param = mxGetPr(prhs[2]);
    double alpha = material_param[0];
    double beta = material_param[1];
    double bulk = material_param[2];
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,Sigma,U_h) private(gradphi,GradUh,ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,alpha,beta,bulk)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        double traceE[NumQuadPoints];
        
        double F[NumQuadPoints][dim][dim];
        double P_Uh[dim][dim];
        double invFT[NumQuadPoints][dim][dim];
        double detF[NumQuadPoints];
        double logdetF[NumQuadPoints];
        double pow2detF[NumQuadPoints];
        double pow23detF[NumQuadPoints];
        double C[NumQuadPoints][dim][dim];
        double I_C[NumQuadPoints];
        
        q = 0;
        
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
        
        detF[q] = MatrixDeterminant(dim, F[q]);
        MatrixInvT(dim, F[q], invFT[q] );
        MatrixProductAlphaT1(dim, 1.0, F[q], F[q], C[q] );
        logdetF[q] = log( detF[q] );
        pow23detF[q] = pow(detF[q], -2.0 / 3.0);
        pow2detF[q] = pow(detF[q], 2.0);
        I_C[q] = Trace(dim, C[q]);
        
        double P1[dim][dim];
        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                P_Uh[d1][d2] =  2.0 * ( alpha + 2.0 * beta * ( pow23detF[q] * I_C[q] - 3.0 ) )
                         * pow23detF[q] * ( F[q][d1][d2] - 1.0 / 3.0 * I_C[q]  * invFT[q][d1][d2] )
                     + 1.0 / 2.0 * bulk * ( pow2detF[q] - detF[q] + logdetF[q] ) * invFT[q][d1][d2];
            }
        }
        
        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                P[ie+(d1+d2*dim)*noe] =  P_Uh[d1][d2] ;
            }
        }
        
        double Sigma_tmp[dim][dim];
        /* Sigma = 1 / det(F) * P * F^T */
        MatrixProductAlphaT2(dim,  1.0 / detF[q], P_Uh, F[q], Sigma_tmp );
        
        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
        {
            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
            {
                Sigma[ie+(d1+d2*dim)*noe] =  Sigma_tmp[d1][d2] ;
            }
        }
    }
    
}
/*************************************************************************/