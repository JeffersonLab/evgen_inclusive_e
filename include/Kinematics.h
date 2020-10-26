#ifndef EVGEN_INCLUSIVE_E_KINEMATICS_H 
#define EVGEN_INCLUSIVE_E_KINEMATICS_H 

// kinematics namespace for common variables like x, Q2, W, etc 
#include <cstdlib>
#include <iostream> 
#include <cmath> 

#include "constants.h"

namespace Kinematics { 
   double GetQ2(double Es,double Ep,double th); 
   double GetW(double Es,double Ep,double th,double M=proton_mass);
   double GetXbj(double Es,double Ep,double th,double M=proton_mass);
   double GetEpsilon(double Es,double Ep,double th);
}

#endif 
