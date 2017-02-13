#ifndef __CHRISTY_BOSTED_INELASTIC_QE_HH
#define __CHRISTY_BOSTED_INELASTIC_QE_HH

//###########################################################################################################
// A direct Translation of Peter Bosted 2009 fortran code to C code
// inherit part of the code from Qweak and Moller(remoll), mainly contributed by Seamus and Ciprian  
// modified to be an isolated C code, can do protton, neutron, A>=2 nucleus
// can reproduce 2007 fortran and 2009 fortran code output
// by Yuxiang Zhao (yxzhao@jlab.org) 
//          Jan 2017
//
//
// F1F2IN09 -- returns the inelastic F1 and F2 structure functions for a nucleus (not a nucleon) with Z and A numbers
// christy507 -- returns F1, R, sigmaT, sigmaL (last two calculated by resmod507_v2)
// resmodd -- returns F1 for average of free proton and neutron
// resmod507_v2 -- calculates sigma T or L
// MEC2009 -- fit to low q2 dip region (purely empirical)
// fitemc -- take into account the EMC effect
//
// based on:
//    http://arxiv.org/pdf/1203.2262v2.pdf --> Empirical fit to e-nucleus scattering.
//    http://arxiv.org/pdf/0712.3731v4.pdf --> e-proton cross sections and resonance fit
//    http://arxiv.org/pdf/0711.0159v4.pdf --> e-D and e-n fits
// data can be found at:
//   http://arxiv.org/pdf/1202.1457v1.pdf -> thesis of V. Mamyan (tables at the back)
//   http://faculty.virginia.edu/qes-archive/QES-data.php --> for QE data
//#######################################################################################################


void christy507(double w2, double q2, double &F1, double &R, double &sigT, double &sigL);
int resmodd(double w2, double q2, double xval[50], double &sig);
double resmod507_v2(double sf, double w2, double q2,double xval[50]);
double MEC2009(double q2, double w2, int A);
double fitemc(double X, int A);

//two additional functions needed for deuteron calculations
double mec(double z, double a, double q2, double w2, double &f1, double &f2, double xval[50]);
double pind(double w2, double q2, double &F1, double &R, double &sigt, double &sigl); //return 0, meaningless
double resd(double q2, double w2, double xval[50]);   //return F1d, per nucleon

int F1F2IN09(int Z, int IA, double qsq, double wsq, double &F1, double &F2, double &Rc);
//returns 0 on success
//returns -1 for failure of resmodd
//returns -2 for A < 3
//returns 1 for nan/inf F1/F2
void F1F2QE09(int Z, int IA, double QSQ, double wsq, double &F1, double &F2);



#endif//__CHRISTY_BOSTED_INELASTIC_QE_HH








