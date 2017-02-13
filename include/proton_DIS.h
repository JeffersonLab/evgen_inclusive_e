#ifndef PROTON_DIS_H
#define PROTON_DIS_H

#include "constants.h"
#include "boost/lexical_cast.hpp"
#include "LHAPDF/LHAPDF.h"
//#include <iostream>
//#include <fstream>
//#include <vector>
#include "TMath.h"
#include "math.h"
#include "TString.h"
using namespace LHAPDF;
using namespace std;

//proton structure functions
double calculate_proton_g1gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g5gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g3gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g5z(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g3z(PDF* pol_pdf, double x, double Q2);

double calculate_proton_F2g(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F1g(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F1gz(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F3gz(PDF* unpol_pdf, double x, double Q2);

//---proton asymmetries
//---PVDIS asymmetry using unpol electron + longitudinally polarized proton
double calculate_proton_AL(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_proton_AL_g1gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_proton_AL_g5gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);

//---PVDIS asymmetry using long. pol. electron + unpolarized proton
double calculate_proton_Abeam(PDF* unpol_pdf, double x, double Q2, double y);


#endif //PROTON_DIS_H
