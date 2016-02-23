


#define nx_inlet  35
#define nx_outlet 35

#define ny_abs 0

#define X_out 570  /**** X_out+nx_in+nx_out ****/  
#define Y_out 60   /**** Y_out+ny_abs ****/
#define Z_out 60     

#define X_m 574   /**** X_out+4 ****/
#define Y_m 64    /**** Y_out+4 ****/
#define Z_m 64      /**** Z_out+4 ****/

#define nx X_out+1    /**** nx+1 ****/
#define ny Y_out+1    /**** ny+1 ****/
#define nz Z_out+1    /**** nz+1 ****/

#define deltaXI 0.0004860290116159057    
#define deltaET 0.008333333333333     //---- 1/(2*(Y_out-ny_abs)) because already normalized ----//
#define deltaZT 0.0003471635797256469

#define nxx nx+1
#define nyy ny+1
#define nzz nz+1

#define nxxx nxx+1
#define nyyy nyy+1
#define nzzz nzz+1 

