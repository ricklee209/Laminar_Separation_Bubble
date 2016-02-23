



#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "Resolution.h"

extern int X_np;

void X_boundary_condition
// ============================================================================ //
(
int myid,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
)
// ============================================================================ //
{
	
#include "ijk.h"
#include "prm.h"


#include "MPI_prm.h"
#include "Mpi_division.h"

	double rho,U,V,W,VV,P,C,T,h,H;
	double tmp0, tmp1, tmp2, temp;

	double XX, YY, eta, Jsum, Jd;
	double p[6];

	p[0] =  0.000231158723443;
	p[1] = -0.001079099290045;
	p[2] = -0.008094664737378;
	p[3] =  0.012993463105787;
	p[4] =  0.326610931292114;
	p[5] =  0.000622120502537;


				
	if (myid == 0) {


		Jsum = Jd = 0.0;

		for (j = 2; j <= ny; j++) Jsum = Jsum + 1./J[3][j][3];


		istart=3;

		for (j = 2; j < nyy; j++) {

			Jd = Jd + 1./J[3][j][3];

			for (k = 2; k < nzz; k++) {

				XX = 1.5*high;
				
				
				eta = (Jd/Jsum)*high/2.0*sqrt( U0/(0.0000185/1.1842*XX) );

				if(eta > 5.0) eta = 5.0;
				

				YY = p[0]*pow(eta,5)\
					+p[1]*pow(eta,4)\
					+p[2]*pow(eta,3)\
					+p[3]*pow(eta,2)\
					+p[4]*pow(eta,1)\
					+p[5]*pow(eta,0);

				V = YY*U0;

				rho = 1.1842;
				U = 0.0;
				W = 0.0;
				P = 101300.0;


				U1_[istart-1][j][k] = rho/J[istart-1][j][k];
				U2_[istart-1][j][k] = rho*U/J[istart-1][j][k];
				U3_[istart-1][j][k] = rho*V/J[istart-1][j][k];
				U4_[istart-1][j][k] = rho*W/J[istart-1][j][k];
				U5_[istart-1][j][k] = (P/(K-1)+0.5*(U*U+V*V+W*W))/J[istart-1][j][k];

			}
		}

	}


	if (myid == np-1) {

		iend = gend[myid];

		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {

				U1_[iend+1][j][k] = U1_[iend][j][k];
				U2_[iend+1][j][k] = U2_[iend][j][k];
				U3_[iend+1][j][k] = U3_[iend][j][k];
				U4_[iend+1][j][k] = U4_[iend][j][k];
				U5_[iend+1][j][k] = U5_[iend][j][k];

			}
		}

	}

}
