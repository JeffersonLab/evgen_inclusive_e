#ifndef norad_XS_H
#define norad_XS_H

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
#include "norad_XS.h"
#include "eInclusiveCrossSection.h"
#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF;
using namespace std;

class noradXS:public eInclusiveCrossSection {
    private: 

    public:

    double GetBornXS(){return calculate_fixed_target_xs( fEs,  fZ,  fA,  fTh,  fEp,  fpdf);}

};
#endif
