#ifndef FIXED_TARGET_XS_H
#define FIXED_TARGET_XS_H

// #include "constants.h"
// #include "boost/lexical_cast.hpp"
// #include "LHAPDF/LHAPDF.h"

#include "constants.h"
#include "boost/lexical_cast.hpp"
#include "LHAPDF/LHAPDF.h"
#include "math.h"

// #include "TMath.h"
// #include "TString.h"

#include "proton_DIS.h"
#include "neutron_DIS.h"
#include "christy_bosted_inelastic_QE.h"
#include "F1F221Wrapper.h"

using namespace LHAPDF;
using namespace std;

double calculate_fixed_target_xs(double E, int Z, int A, double theta, double Ep, PDF* unpol_pdf, int Fit_model);

#endif

