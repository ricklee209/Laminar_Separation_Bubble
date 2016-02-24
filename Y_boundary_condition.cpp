




#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "Resolution.h"

extern int X_np;


#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


void Y_boundary_condition
(
 // ============================================================================ //
 int myid,

 double deltaT,

 double heat_flux,

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

 double (*etdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

 double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]

 // ============================================================================ //
 )

{
	
#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"
	
	double rho,U,V,W,VV,P,C,T,h,H;
	double temp;

//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 0;            			  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+3;				  ////
//// ============================================ ////

//#pragma omp parallel for private(k,rho,U,V,W,VV,P,temp,T)
	for (i = istart; i <= iend; i++) {
		for (k = 2; k <= nz; k++) {
			
			/*
			U1_[i][1][k] = U1_[i][2][k]*J[i][2][k]/J[i][1][k];
			U2_[i][1][k] = -U2_[i][2][k]*J[i][2][k]/J[i][1][k];
			U3_[i][1][k] = -U3_[i][2][k]*J[i][2][k]/J[i][1][k];
			U4_[i][1][k] = -U4_[i][2][k]*J[i][2][k]/J[i][1][k];
			U5_[i][1][k] = U5_[i][2][k]*J[i][2][k]/J[i][1][k];

			U1_[i][0][k] = U1_[i][3][k]*J[i][3][k]/J[i][0][k];
			U2_[i][0][k] = -U2_[i][3][k]*J[i][3][k]/J[i][0][k];
			U3_[i][0][k] = -U3_[i][3][k]*J[i][3][k]/J[i][0][k];
			U4_[i][0][k] = -U4_[i][3][k]*J[i][3][k]/J[i][0][k];
			U5_[i][0][k] = U5_[i][3][k]*J[i][3][k]/J[i][0][k];


			U1_[i][nyy][k] = U1_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U2_[i][nyy][k] = -U2_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U3_[i][nyy][k] = -U3_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U4_[i][nyy][k] = -U4_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U5_[i][nyy][k] = U5_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];

			U1_[i][nyyy][k] = U1_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U2_[i][nyyy][k] = -U2_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U3_[i][nyyy][k] = -U3_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U4_[i][nyyy][k] = -U4_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U5_[i][nyyy][k] = U5_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			*/

			rho = U1_[i][2][k]*J[i][2][k];
			U = U2_[i][2][k]/U1_[i][2][k];
			V = U3_[i][2][k]/U1_[i][2][k];
			W = U4_[i][2][k]/U1_[i][2][k];     
			VV = U*U+V*V+W*W;
			P = (U5_[i][2][k]*J[i][2][k]-0.5*rho*VV)*(K-1);
			temp = P/rho/R;

			T = 2*298.0592-temp;

			rho = P/R/T;

			U1_[i][1][k] = P/R/T/J[i][1][k];
			U2_[i][1][k] = -P/R/T*U/J[i][1][k];
			U3_[i][1][k] = -P/R/T*V/J[i][1][k];
			U4_[i][1][k] = -P/R/T*W/J[i][1][k];
			U5_[i][1][k] = (rho*R*T/(K-1)+0.5*P/R/T*VV)/J[i][1][k];




			rho = U1_[i][ny][k]/J[i][ny][k];
			U = U2_[i][ny][k]/U1_[i][ny][k];

			temp = ( (i+gstart[myid])*deltaXI-1.5*high )/(0.12*high);

			V = 0.7*U0*exp(-temp*temp);

			//if(k == 2) printf("%d\t%f\t%f\n",myid,V,(i+gstart[myid])*1.0);

			W = U3[i][ny][k]/U1_[i][ny][k];   
			VV = U*U+V*V+W*W;
			P = (U5_[i][2][k]*J[i][2][k]-0.5*rho*VV)*(K-1);


			U1_[i][nyy][k] = U1_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U2_[i][nyy][k] = U2_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U3_[i][nyy][k] = rho*V/J[i][ny][k];
			U4_[i][nyy][k] = U2_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U5_[i][nyy][k] = ( P/(K-1.0)+0.5*rho*VV )/J[i][nyy][k];

			U1_[i][nyyy][k] = U1_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U2_[i][nyyy][k] = -U2_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U3_[i][nyyy][k] = -U3_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U4_[i][nyyy][k] = -U4_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U5_[i][nyyy][k] = U5_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			
		}
	}
#pragma omp barrier



}
