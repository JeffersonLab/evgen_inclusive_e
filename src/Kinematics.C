#include "Kinematics.h"
//______________________________________________________________________________
namespace Kinematics {
   //________________________________________________________________________
   double GetQ2(double Es,double Ep,double th){
      double thr  = th*deg_to_rad;
      double SIN  = sin(thr/2.);
      double SIN2 = SIN*SIN;
      double Q2   = 4.*Es*Ep*SIN2;
      return Q2;
   }
   //________________________________________________________________________
   double GetW(double Es,double Ep,double th){
      double Nu = Es-Ep;
      double Q2 = GetQ2(Es,Ep,th);
      double W2 = proton_mass*proton_mass + 2.*proton_mass*Nu - Q2;
      return sqrt(W2);
   }
}
