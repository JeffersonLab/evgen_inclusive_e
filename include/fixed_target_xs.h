#ifndef FIXED_TARGET_XS_H
#define FIXED_TARGET_XS_H

#include "constants.h"
#include "boost/lexical_cast.hpp"
#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF;
using namespace std;



double calculate_fixed_target_xs(double E, int Z, int A, double theta, double Ep, PDF* unpol_pdf);


#endif

