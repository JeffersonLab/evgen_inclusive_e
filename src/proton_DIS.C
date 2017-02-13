

#include "constants.h"
#include "boost/lexical_cast.hpp"
#include "LHAPDF/LHAPDF.h"
//#include <iostream>
//#include <fstream>
//#include <vector>
#include "TMath.h"
#include "math.h"
#include "TString.h"

#include "proton_DIS.h"

using namespace LHAPDF;
using namespace std;

/*
//---functions defined in proton_DIS.h



//proton structure functions
//---pol structure functions
double calculate_proton_g1gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g5gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g3gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g5z(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g3z(PDF* pol_pdf, double x, double Q2);



//---unpol structure functions
double calculate_proton_F2g(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F1g(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F1gz(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F3gz(PDF* unpol_pdf, double x, double Q2);


//---proton asymmetries
//---PVDIS asymmetry with unpolarized electron and longitudinally polarized proton
double calculate_proton_Abeam(PDF* unpol_pdf, double x, double Q2, double y);


double calculate_proton_AL(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_proton_AL_g1gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_proton_AL_g5gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
*/


//----------------------------------------------------------------------
// several functions to calculate structure functions, asymmetries
// only use u, d, s       No c, b, t





//#################################################################
double calculate_proton_g1gz(PDF* pol_pdf, double x, double Q2){
	// g1gz=\sum e_q (gv_q) * (delt_q + delta_qbar)
	
	// delta_u
	double delta_u=pol_pdf->xfxQ2(2, x, Q2); 
	delta_u=delta_u/x;

	//delta_ubar
	double delta_ubar=pol_pdf->xfxQ2(-2, x, Q2);
	delta_ubar=delta_ubar/x;

	//delta_d
	double delta_d=pol_pdf->xfxQ2(1, x, Q2); 
	delta_d=delta_d/x;

	//delta_dbar
	double delta_dbar=pol_pdf->xfxQ2(-1, x, Q2);
	delta_dbar=delta_dbar/x;
	
	//delta_s
	double delta_s=pol_pdf->xfxQ2(3, x, Q2); 
	delta_s=delta_s/x;
	
	//delta_sbar
	double delta_sbar=pol_pdf->xfxQ2(-3, x, Q2); 
	delta_sbar=delta_sbar/x;

	double g1gz= Q_u * gv_u * (delta_u + delta_ubar)
	           + Q_d * gv_d * (delta_d + delta_dbar)
		   + Q_s * gv_s * (delta_s + delta_sbar);
/*	
	cout<<"debug: delta_u= "<<delta_u<<endl;
	cout<<"debug: delta_d= "<<delta_d<<endl;
	cout<<"debug: delta_s= "<<delta_s<<endl;
	cout<<"debug: delta_ubar= "<<delta_ubar<<endl;
	cout<<"debug: delta_dbar= "<<delta_dbar<<endl;
	cout<<"debug: delta_sbar= "<<delta_sbar<<endl;
*/
	return g1gz;

}


//##############################################################
double calculate_proton_g5gz(PDF* pol_pdf, double x, double Q2){
	// g5gz=\sum e_q (gA_q) * (delt_q - delta_qbar)
	
	// delta_u
	double delta_u=pol_pdf->xfxQ2(2, x, Q2); 
	delta_u=delta_u/x;

	//delta_ubar
	double delta_ubar=pol_pdf->xfxQ2(-2, x, Q2);
	delta_ubar=delta_ubar/x;

	//delta_d
	double delta_d=pol_pdf->xfxQ2(1, x, Q2); 
	delta_d=delta_d/x;

	//delta_dbar
	double delta_dbar=pol_pdf->xfxQ2(-1, x, Q2);
	delta_dbar=delta_dbar/x;
	
	//delta_s
	double delta_s=pol_pdf->xfxQ2(3, x, Q2); 
	delta_s=delta_s/x;
	
	//delta_sbar
	double delta_sbar=pol_pdf->xfxQ2(-3, x, Q2); 
	delta_sbar=delta_sbar/x;

	double g5gz=  Q_u * gA_u * (delta_u - delta_ubar)
	            + Q_d * gA_d * (delta_d - delta_dbar)
		    + Q_s * gA_s * (delta_s - delta_sbar);

	return g5gz;
}

//###############################################################
double calculate_proton_F2g(PDF* unpol_pdf, double x, double Q2){
	
	//u
	double u=unpol_pdf->xfxQ2(2, x, Q2);
	u=u/x;

	//ubar
	double ubar=unpol_pdf->xfxQ2(-2, x, Q2);
	ubar=ubar/x;
	
	//d
	double d=unpol_pdf->xfxQ2(1, x, Q2);
	d=d/x;

	//dbar
	double dbar=unpol_pdf->xfxQ2(-1, x, Q2);
	dbar=dbar/x;
	
	//s
	double s=unpol_pdf->xfxQ2(3, x, Q2);
	s=s/x;

	//sbar
	double sbar=unpol_pdf->xfxQ2(-3, x, Q2);
	sbar=sbar/x;

	double F2g=  Q_u * Q_u * (u + ubar) * x
		   + Q_d * Q_d * (d + dbar) * x
		   + Q_s * Q_s * (s + sbar) * x;
	
	return F2g;
}


//###############################################################
double calculate_proton_F1g(PDF* unpol_pdf, double x, double Q2){
	// F1g = 1/2 \sum e_q^2 (q + qbar)
	
	//u
	double u=unpol_pdf->xfxQ2(2, x, Q2);
	u=u/x;

	//ubar
	double ubar=unpol_pdf->xfxQ2(-2, x, Q2);
	ubar=ubar/x;
	
	//d
	double d=unpol_pdf->xfxQ2(1, x, Q2);
	d=d/x;

	//dbar
	double dbar=unpol_pdf->xfxQ2(-1, x, Q2);
	dbar=dbar/x;
	
	//s
	double s=unpol_pdf->xfxQ2(3, x, Q2);
	s=s/x;

	//sbar
	double sbar=unpol_pdf->xfxQ2(-3, x, Q2);
	sbar=sbar/x;

	double F2g=  Q_u * Q_u * (u + ubar) * x
		   + Q_d * Q_d * (d + dbar) * x
		   + Q_s * Q_s * (s + sbar) * x;
	
	//---to calculate F1g, one need R corrections
	// calculate F2g first then do F1g calculation
	//hep-ex/9808028 (1998)
	double a1=0.0485, a2=0.5470, a3=2.0621, a4=-0.3804, a5=0.5090, a6=-0.0285;
	double b1=0.0481, b2=0.6114, b3=-0.3509,b4=-0.4611, b5=0.7172, b6=-0.0317;
	double c1=0.0577, c2=0.4644, c3=1.8288, c4=12.3708, c5=-43.1043,c6=41.7415;

	double theta=1.0 + 12.0*(Q2/(Q2+1.0))*(0.125*0.125/(0.125*0.125+x*x));

	double Q2_thr=c4*x + c5*x*x + c6*x*x*x;

	double R_const=theta/(TMath::Log(Q2/0.04));

	double R_a = a1*R_const + a2/pow( pow(Q2, 4) + pow(a3,4), 1.0/4.0)  * ( 1.0+a4*x+a5*x*x )* pow(x, a6);
	double R_b = b1*R_const + (b2/Q2 + b3/(Q2*Q2 + 0.3*0.3))*(1.0+b4*x+b5*x*x)*pow(x,b6);
	double R_c = c1*R_const + c2*pow( pow(Q2-Q2_thr, 2) + c3*c3    , -0.5);

	double R_average=(R_a + R_b + R_c)/3.0;

//	gamma2=4.*AMP**2.*xbj**2./qmu2

	double MP=0.93827;
	double gamma2=4.0*MP*MP*x*x/Q2;
	
//	f1p=f2p*(1+gamma2)/(2*xbj*(1+r1998(xbj,qmu2)))

	double F1g=F2g*(1.0+gamma2)/(2.0*x*(1.0+R_average));
/*	
	cout<<"debug: u= "<<u<<"	d= "<<d<<"	s= "<<s<<endl;
	cout<<"debug: ubar= "<<ubar<<"	dbar= "<<dbar<<"	sbar= "<<sbar<<endl;
	cout<<"debug: R= "<<R_average<<endl;
	cout<<"debug: F1g= "<<F1g<<"    F2= "<<F2g<<endl;
	cout<<"debug: gamma2= "<<gamma2<<" F1g/F2g= "<<(1.0+gamma2)/(2.0*x*(1.0+R_average))<<endl;
*/	return F1g;
}


//###############################################################
double calculate_proton_g3gz(PDF* pol_pdf, double x, double Q2){
	double g5gz=calculate_proton_g5gz(pol_pdf, x, Q2);
	double g3gz=2.0 * x * g5gz;
	return g3gz;
}


//#############################################################
double calculate_proton_g5z(PDF* pol_pdf, double x, double Q2){
	// g5z=\sum (gv_q) (gA_q) * (delt_q - delta_qbar)
	
	// delta_u
	double delta_u=pol_pdf->xfxQ2(2, x, Q2); 
	delta_u=delta_u/x;

	//delta_ubar
	double delta_ubar=pol_pdf->xfxQ2(-2, x, Q2);
	delta_ubar=delta_ubar/x;

	//delta_d
	double delta_d=pol_pdf->xfxQ2(1, x, Q2); 
	delta_d=delta_d/x;

	//delta_dbar
	double delta_dbar=pol_pdf->xfxQ2(-1, x, Q2);
	delta_dbar=delta_dbar/x;
	
	//delta_s
	double delta_s=pol_pdf->xfxQ2(3, x, Q2); 
	delta_s=delta_s/x;
	
	//delta_sbar
	double delta_sbar=pol_pdf->xfxQ2(-3, x, Q2); 
	delta_sbar=delta_sbar/x;

	double g5z= gv_u * gA_u * (delta_u - delta_ubar)
	          + gv_d * gA_d * (delta_d - delta_dbar)
		  + gv_s * gA_s * (delta_s - delta_sbar);
	
	return g5z;
}


//#############################################################
double calculate_proton_g3z(PDF* pol_pdf, double x, double Q2){
	double g5z=calculate_proton_g5z(pol_pdf, x, Q2);
	double g3z=2.0 * x * g5z;
	return g3z;
}


//#######################################################################################
double calculate_proton_AL(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y){

	double prob_const = GF*Q2/2.0/sqrt(2.0)/PI/alpha;
	double Y_m = 2.0*y - y*y;
	double Y_p = y*y -2.0*y + 2;

	double F1g = calculate_proton_F1g(unpol_pdf, x, Q2);
	double g1gz = calculate_proton_g1gz(pol_pdf, x, Q2);
	double g5gz = calculate_proton_g5gz(pol_pdf, x, Q2);

	double AL= prob_const * ( gv_e*g5gz/F1g + gA_e*Y_m/Y_p*g1gz/F1g );

	return AL;
}

//############################################################################################
double calculate_proton_AL_g1gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y){

	double prob_const = GF*Q2/2.0/sqrt(2.0)/PI/alpha;
	double Y_m = 2.0*y - y*y;
	double Y_p = y*y -2.0*y + 2;

	double F1g = calculate_proton_F1g(unpol_pdf, x, Q2);
	double g1gz = calculate_proton_g1gz(pol_pdf, x, Q2);
	double g5gz = calculate_proton_g5gz(pol_pdf, x, Q2);

	double AL_g1gz= prob_const * (  gA_e*Y_m/Y_p*g1gz/F1g );

	return AL_g1gz;
}

//###########################################################################################
double calculate_proton_AL_g5gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y){

	double prob_const = GF*Q2/2.0/sqrt(2.0)/PI/alpha;
	double Y_m = 2.0*y - y*y;
	double Y_p = y*y -2.0*y + 2;

	double F1g = calculate_proton_F1g(unpol_pdf, x, Q2);
	double g1gz = calculate_proton_g1gz(pol_pdf, x, Q2);
	double g5gz = calculate_proton_g5gz(pol_pdf, x, Q2);

	double AL_g5gz= prob_const * ( gv_e*g5gz/F1g );

	return AL_g5gz;
}


//####################################################################
double calculate_proton_F1gz(PDF* unpol_pdf, double x, double Q2){
	
	//u
	double u=unpol_pdf->xfxQ2(2, x, Q2);
	u=u/x;

	//ubar
	double ubar=unpol_pdf->xfxQ2(-2, x, Q2);
	ubar=ubar/x;
	
	//d
	double d=unpol_pdf->xfxQ2(1, x, Q2);
	d=d/x;

	//dbar
	double dbar=unpol_pdf->xfxQ2(-1, x, Q2);
	dbar=dbar/x;
	
	//s
	double s=unpol_pdf->xfxQ2(3, x, Q2);
	s=s/x;

	//sbar
	double sbar=unpol_pdf->xfxQ2(-3, x, Q2);
	sbar=sbar/x;

	double F1gz=  Q_u * gv_u * (u+ubar)
		    + Q_d * gv_d * (d+dbar)
		    + Q_s * gv_s * (s+sbar);
	
	return F1gz;
}


//###################################################################
double calculate_proton_F3gz(PDF* unpol_pdf, double x, double Q2){
	
	//u
	double u=unpol_pdf->xfxQ2(2, x, Q2);
	u=u/x;

	//ubar
	double ubar=unpol_pdf->xfxQ2(-2, x, Q2);
	ubar=ubar/x;
	
	//d
	double d=unpol_pdf->xfxQ2(1, x, Q2);
	d=d/x;

	//dbar
	double dbar=unpol_pdf->xfxQ2(-1, x, Q2);
	dbar=dbar/x;
	
	//s
	double s=unpol_pdf->xfxQ2(3, x, Q2);
	s=s/x;

	//sbar
	double sbar=unpol_pdf->xfxQ2(-3, x, Q2);
	sbar=sbar/x;

	double F3gz=  Q_u * gA_u * (u - ubar)
		    + Q_d * gA_d * (d - dbar)
		    + Q_s * gA_s * (s - sbar);
	
	F3gz = F3gz*2.0;

	return F3gz;
}

//#######################################################################
double calculate_proton_Abeam(PDF* unpol_pdf, double x, double Q2, double y){
	double F1g = calculate_proton_F1g(unpol_pdf, x, Q2);
	double F1gz = calculate_proton_F1gz(unpol_pdf, x, Q2);
	double F3gz = calculate_proton_F3gz(unpol_pdf, x, Q2);
	double Abeam= GF*Q2/2/sqrt(2)/PI/alpha * ( gA_e * F1gz/F1g + gv_e * (2*y-y*y)/2.0/(y*y-2*y+2) * F3gz/F1g  );
	return Abeam;
}
