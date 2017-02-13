

#include "constants.h"
#include "boost/lexical_cast.hpp"
#include "LHAPDF/LHAPDF.h"
#include "TMath.h"
#include "math.h"
#include "TString.h"

#include "proton_DIS.h"
#include "neutron_DIS.h"
#include "christy_bosted_inelastic_QE.h"
#include "fixed_target_xs.h"


using namespace LHAPDF;
using namespace std;

double calculate_fixed_target_xs(double E, int Z, int A, double theta, double Ep, PDF* unpol_pdf){
	// E : electron beam energy, in unit of GeV
	// Z : target z
	// A : target A, N=A-Z
	// theta : outgoing electron angle, in unit of degree, range: [0,180] degree
	// Ep : E' of the outgoing electron, in unit of GeV
	//
	//return nucleus, NOT nucleon, xs value in unit of mub/GeV-sr
	
	double Nu=E-Ep;
	double Q2=4.0*E*Ep*sin(theta*deg_to_rad/2.0)*sin(theta*deg_to_rad/2.0);
	double x=Q2/2.0/proton_mass/Nu;
	double y=Nu/E;
	double W2=proton_mass * proton_mass + 2 * proton_mass * Nu - Q2;
	
	double F1, F2, xs, r;
	if(W2<3.0*3.0){  //use PB model
		F1F2IN09(Z, A, Q2, W2, F1, F2, r);
	}else if(W2>=3.0*3.0){ //use pdf sets
		double F1p=calculate_proton_F1g(unpol_pdf,  x,  Q2);  //proton 
		double F2p=calculate_proton_F2g(unpol_pdf,  x,  Q2);

		double F1n=calculate_neutron_F1g( unpol_pdf,  x, Q2); //neutron
		double F2n=calculate_neutron_F2g( unpol_pdf,  x,  Q2);

		F1=Z*F1p + (A-Z)*F1n;  //summing over proton and neutron
		F2=Z*F2p + (A-Z)*F2n;
	}

	//already get F1 and F2 for the nucleus
	xs=(2./137.*Ep/Q2*cos(theta*deg_to_rad/2))*(2./137.*Ep/Q2*cos(theta*deg_to_rad/2));     // mott
	xs=xs*(2/proton_mass*F1*tan(abs(theta*deg_to_rad)/2)*tan(abs(theta*deg_to_rad)/2) + F2/Nu);
	xs=xs*389.379;   //  mub/GeV-sr

	return xs;
}




		
		
		



