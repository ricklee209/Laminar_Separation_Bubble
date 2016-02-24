



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int X_np;

void RK3NODT
(
// ============================================================================ //
int myid,

int RK,

double deltaT,

double *er1,
double *er2,
double *er3,
double *er4,
double *er5,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*U1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*U1q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*xidx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFx1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*inFy1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*inFz1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*vF2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*vF3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*vF4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*vF5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 


/**** MR = RHS ****/
double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
/**** MR = RHS-end ****/


// ============================================================================ //
)

{

	
#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"

	double tmp;

	double Rp1,Rp2,Rp3,Rp4,Rp5,
		   Rf1,Rf2,Rf3,Rf4,Rf5,
		   Rk1,Rk2,Rk3,Rk4,Rk5;

	double c1[3] = { 1.0, 0.75, 1.0/3.0};
	double c2[3] = { 0.0, 0.25, 2.0/3.0};
	double c3[3] = { 1.0, 0.25, 2.0/3.0};





	
// ------------------------------- Absorbing boundary condition ------------------------------- //
// -------------------------------------------------------------------------------------------- //

	double C_plan;

	double U_in_0 = 1.15;
	
	double Sigma_in_0 = 0.035;

	double U_out_0 = 1.15;
	
	double Sigma_out_0 = 1.25;
	
	double U_in_1, U_out_1;
	
	double Sigma_in, Sigma_out;

	double rho,U,V,W,VV,P,C,T,H;
	double u,v,w,beta;

	if( (gstart[myid]) <= nx_inlet-1 ) {
	
//// ============================================ ////
		istart =  3;	
//// ============================================ ////
		if ( (gend0[myid]) <= nx_inlet-1 )
			iend = gend[myid];
		else
			iend = nx_inlet-gstart[myid]+2;
//// ============================================ ////


	for (i = istart ; i <= iend; i++) {
			
#pragma omp parallel for private(C,C_plan,beta,k,rho,u,v,w,VV,P,U_in_1,Sigma_in)

			for (j = 2; j <= ny; j++) {
				for (k = 2; k <= nz; k++) {

					rho = U1_[i][j][k]*J[i][j][k];
					u = U2_[i][j][k]/U1_[i][j][k];
					v = U3_[i][j][k]/U1_[i][j][k];
					w = U4_[i][j][k]/U1_[i][j][k];     
					VV = u*u+v*v+w*w;
					P = (U5_[i][j][k]*J[i][j][k]-0.5*rho*VV)*(K-1);
					C = K*P/rho;

					/* preconditioning */
					beta = max(VV/C,e);

					C_plan = 0.5*sqrt(u*u*(beta-1)*(beta-1)+4*beta*C);
					
					U_in_1 = pow( ((nx_inlet-(i-2-gstart[myid]-0.5)*1.0)/nx_inlet), 3.0 )* U_in_0*C_plan;
					
					Sigma_in = Sigma_in_0*U_in_1/high;

					
					U1q[i][j][k] = -U_in_1*( U1_[i][j][k]*J[i][j][k]-U1_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_in*(rho-1.1842);
					U2q[i][j][k] = -U_in_1*( U2_[i][j][k]*J[i][j][k]-U2_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_in*(rho*u-1.1842*U0);
					U3q[i][j][k] = -U_in_1*( U3_[i][j][k]*J[i][j][k]-U3_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_in*(rho*v);
					U4q[i][j][k] = -U_in_1*( U4_[i][j][k]*J[i][j][k]-U4_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_in*(rho*w);
					U5q[i][j][k] = -U_in_1*( U5_[i][j][k]*J[i][j][k]-U5_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_in*(U5_[i][j][k]*J[i][j][k]-253250.0-0.5*1.1842*U0*U0);


					//if(j ==2 && k == 2) printf("%f\t%f\t%f\n",U_in_1, Sigma_in,high);
					
				}
			}
		}

#pragma omp barrier

	}


	
	if ( gend0[myid] >= (X_out-nx_outlet) ) {

//// ============================================ ////

		if ( (gstart[myid]) >= (X_out-nx_outlet))
			istart =  3;	
		else 
			istart = (X_out-nx_outlet)+3-gstart[myid];

//// ============================================ ////
		iend = gend[myid];			    		  ////
//// ============================================ ////


		for (i = istart ; i <= iend; i++) {
			
#pragma omp parallel for private(C,C_plan,beta,k,rho,u,v,w,VV,P,U_out_1,Sigma_out)

			for (j = 2; j <= ny; j++) {
				for (k = 2; k <= nz; k++) {

					rho = U1_[i][j][k]*J[i][j][k];
					u = U2_[i][j][k]/U1_[i][j][k];
					v = U3_[i][j][k]/U1_[i][j][k];
					w = U4_[i][j][k]/U1_[i][j][k];     
					VV = u*u+v*v+w*w;
					P = (U5_[i][j][k]*J[i][j][k]-0.5*rho*VV)*(K-1);
					C = K*P/rho;

					/* preconditioning */
					beta = max(VV/C,e);

					C_plan = 0.5*sqrt(u*u*(beta-1)*(beta-1)+4*beta*C);
					
					U_out_1 = pow( (((i-2+gstart[myid]-0.5-(X_out-nx_outlet))*1.0)/nx_outlet), 3.0 )* U_out_0*C_plan;

					Sigma_out = Sigma_out_0*U_out_1/high;
					
					
					
					U1q[i][j][k] = -U_out_1*( U1_[i][j][k]*J[i][j][k]-U1_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_out*(rho-1.1842);
					U2q[i][j][k] = -U_out_1*( U2_[i][j][k]*J[i][j][k]-U2_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_out*(rho*u-1.1842*U0);
					U3q[i][j][k] = -U_out_1*( U3_[i][j][k]*J[i][j][k]-U3_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_out*(rho*v);
					U4q[i][j][k] = -U_out_1*( U4_[i][j][k]*J[i][j][k]-U4_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_out*(rho*w);
					U5q[i][j][k] = -U_out_1*( U5_[i][j][k]*J[i][j][k]-U5_[i-1][j][k]*J[i-1][j][k] )/deltaXI-Sigma_out*(U5_[i][j][k]*J[i][j][k]-253250.0-0.5*1.1842*U0*U0);
					
					
					//if(j ==2 && k == 2) printf("%f\t%f\t%f\n",U_out_1, Sigma_out,C_plan);
					
					
				}
			}
		}

#pragma omp barrier

	}

// -------------------------------------------------------------------------------------------- //
// ------------------------------- Absorbing boundary condition ------------------------------- //




	

//// ============================================ ////
			istart = 3; 	            		  ////	
//// ============================================ ////
			iend = gend[myid];		    		  ////
//// ============================================ ////


	for (i = istart ; i <= iend; i++) {
		
#pragma omp parallel for private(\
k,_k,tmp,\
Rf1,Rf2,Rf3,Rf4,Rf5,Rk1,Rk2,Rk3,Rk4,Rk5\
)
		for (j = 2; j < nyy; j++) {
			for (k = 2, _k = 1; k < nzz; k++, _k++) {

				tmp = U2_[i][j][k]/U1_[i][j][k];

				Rf1 = -((inFx1[i][j-1][_k]-inFx1[i-1][j-1][_k])/deltaXI+
					(inFy1[i-1][j][_k]-inFy1[i-1][j-1][_k])/deltaET+
					(inFz1[i-1][j-1][k]-inFz1[i-1][j-1][_k])/deltaZT);
					
				Rf2 = -((inFx2[i][j-1][_k]-inFx2[i-1][j-1][_k])/deltaXI+
					(inFy2[i-1][j][_k]-inFy2[i-1][j-1][_k])/deltaET+
					(inFz2[i-1][j-1][k]-inFz2[i-1][j-1][_k])/deltaZT);
					
					
				Rf3 = -((inFx3[i][j-1][_k]-inFx3[i-1][j-1][_k])/deltaXI+
					(inFy3[i-1][j][_k]-inFy3[i-1][j-1][_k])/deltaET+
					(inFz3[i-1][j-1][k]-inFz3[i-1][j-1][_k])/deltaZT);
					
				Rf4 = -((inFx4[i][j-1][_k]-inFx4[i-1][j-1][_k])/deltaXI+
					(inFy4[i-1][j][_k]-inFy4[i-1][j-1][_k])/deltaET+
					(inFz4[i-1][j-1][k]-inFz4[i-1][j-1][_k])/deltaZT);
					
				Rf5 = -((inFx5[i][j-1][_k]-inFx5[i-1][j-1][_k])/deltaXI+
					(inFy5[i-1][j][_k]-inFy5[i-1][j-1][_k])/deltaET+
					(inFz5[i-1][j-1][k]-inFz5[i-1][j-1][_k])/deltaZT);


				Rk1 = U1q[i][j][k]/J[i][j][k]+Rf1;
				Rk2 = U2q[i][j][k]/J[i][j][k]+Rf2+vF2[i][j][k];
				Rk3 = U3q[i][j][k]/J[i][j][k]+Rf3+vF3[i][j][k];
				Rk4 = U4q[i][j][k]/J[i][j][k]+Rf4+vF4[i][j][k];
				Rk5 = U5q[i][j][k]/J[i][j][k]+Rf5+vF5[i][j][k];


				MR1[i][j][k] = deltaT*Rk1;
				MR2[i][j][k] = deltaT*Rk2;
				MR3[i][j][k] = deltaT*Rk3;
				MR4[i][j][k] = deltaT*Rk4;
				MR5[i][j][k] = deltaT*Rk5;
				

				U1_[i][j][k] = c1[RK-1]*U1[i][j][k] + c2[RK-1]*U1_[i][j][k] + c3[RK-1]*MR1[i][j][k];
				U2_[i][j][k] = c1[RK-1]*U2[i][j][k] + c2[RK-1]*U2_[i][j][k] + c3[RK-1]*MR2[i][j][k];
				U3_[i][j][k] = c1[RK-1]*U3[i][j][k] + c2[RK-1]*U3_[i][j][k] + c3[RK-1]*MR3[i][j][k];
				U4_[i][j][k] = c1[RK-1]*U4[i][j][k] + c2[RK-1]*U4_[i][j][k] + c3[RK-1]*MR4[i][j][k];
				U5_[i][j][k] = c1[RK-1]*U5[i][j][k] + c2[RK-1]*U5_[i][j][k] + c3[RK-1]*MR5[i][j][k];

			}
		}
	}
#pragma omp barrier



}