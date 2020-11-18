#ifndef DIPOLE_ELASTIC_FORM_FACTOR_HH
#define DIPOLE_ELASTIC_FORM_FACTOR_HH

// dipole model of the elastic form factor for the proton 
// fit from GMp(12) preprint  

#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "ElasticFormFactor.h"

class DipoleElasticFF: public ElasticFormFactor { 

   private:

   public: 
      DipoleElasticFF();
      ~DipoleElasticFF();

      double GetGE(double Q2); 
      double GetGM(double Q2); 

}; 

#endif 
