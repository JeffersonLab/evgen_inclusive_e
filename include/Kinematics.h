#ifndef EVGEN_INCLUSIVE_E_KINEMATICS_H 
#define EVGEN_INCLUSIVE_E_KINEMATICS_H 

// kinematics namespace for common variables like x, Q2, W, etc 
#include <cmath> 

#include "constants.h"

namespace Kinematics { 
   double GetQ2(double Es,double Ep,double th); 
   double GetW(double Es,double Ep,double th); 
}

#endif 
