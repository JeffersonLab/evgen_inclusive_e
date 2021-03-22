#ifndef eINCLUSIVECROSSSECTION_H
#define eINCLUSIVECROSSSECTION_H

// eInclusiveCrossSection class
// Notes: - Abstract base class 
//        - Variables: beam energy (Es) and scattered energy (Ep) are in GeV, 
//          scattering angle (th) is in degrees.  A is in g/mol. 
//        - Needed values: Z, A, Es, Ep, th 
//        - The resulting cross section is for inclusive inelastic electron scattering  

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "constants.h"
#include "boost/lexical_cast.hpp"
#include "LHAPDF/LHAPDF.h"
#include "math.h"

// #include "TMath.h"
// #include "TString.h"

#include "proton_DIS.h"
#include "neutron_DIS.h"
#include "christy_bosted_inelastic_QE.h"
#include "fixed_target_xs.h"
#include "eInclusiveCrossSection.h"
#include "LHAPDF/LHAPDF.h"

using namespace LHAPDF;
using namespace std;

class eInclusiveCrossSection{

   protected:
      double fZ,fA;
      double fEs,fEp,fTh;
      PDF *fpdf; 
      double fscale;
      double fmodel;
      void Init();

   public:
      eInclusiveCrossSection();
      ~eInclusiveCrossSection();

      void SetZ(double v)  { fZ  = v; }
      void SetA(double v)  { fA  = v; }
      void SetTh(double v) { fTh = v; }
      void SetEs(double v) { fEs = v; }
      void SetEp(double v) { fEp = v; }
      void Setpdf(PDF* v)  { fpdf= v; }
      void Setmodel(double v)  { fmodel= v; }
      void SetScale(double v)  { fscale= v; }
      double GetZ()  const { return fZ;  }
      double GetA()  const { return fA;  }
      double GetEs() const { return fEs; }
      double GetEp() const { return fEp; }
      double GetTh() const { return fTh; }
      double GetMottXS(double,double);

      virtual double GetBornXS(){return 0;}

};

#endif 
