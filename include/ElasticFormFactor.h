#ifndef ELASTIC_FORM_FACTOR_H
#define ELASTIC_FORM_FACTOR_H

// Elastic form factor class
// Notes: - Abstract base class 
//        - User provides a model deriving from this to determine GE and GM as a function of Q2   

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "constants.h"

class ElasticFormFactor {

   protected:
      void Init();

   public:
      ElasticFormFactor();
      ~ElasticFormFactor();

      virtual double GetGE(double Q2)=0;
      virtual double GetGM(double Q2)=0;

};

#endif 
