#ifndef NEUTRON_DIS_H
#define NEUTRON_DIS_H



#include "constants.h"
#include "boost/lexical_cast.hpp"
#include "LHAPDF/LHAPDF.h"
#include "TMath.h"
#include "math.h"
#include "TString.h"
using namespace LHAPDF;
using namespace std;


//---structure functions
//using isospin to get it from PDF sets
double calculate_neutron_g1gz(PDF* pol_pdf, double x, double Q2);

double calculate_neutron_g5gz(PDF* pol_pdf, double x, double Q2);

double calculate_neutron_F2g(PDF* unpol_pdf, double x, double Q2);

double calculate_neutron_F1g(PDF* unpol_pdf, double x, double Q2);

double calculate_neutron_F1gz(PDF* unpol_pdf, double x, double Q2);

double calculate_neutron_F3gz(PDF* unpol_pdf, double x, double Q2);


//---PVDIS asymmetries

double calculate_neutron_Abeam(PDF* unpol_pdf, double x, double Q2, double y);

double calculate_neutron_AL(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);

double calculate_neutron_AL_g1gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);

double calculate_neutron_AL_g5gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);



#endif   //NEUTRON_DIS_H
