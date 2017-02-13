#ifndef CONSTANTS_YX_H
#define CONSTANTS_YX_H

const double sin2thetaW=0.232;

const double gv_e=-0.5+2.0*sin2thetaW;
const double gA_e=-0.5;

// u
const double gv_u=0.5-4.0/3.0*sin2thetaW;
const double gA_u=0.5;

//d
const double gv_d=-0.5+2.0/3.0*sin2thetaW;
const double gA_d=-0.5;

//s
const double gv_s=-0.5+2.0/3.0*sin2thetaW;
const double gA_s=-0.5;

//charge
const double Q_e=-1;
const double Q_u=2.0/3.0;
const double Q_d=-1.0/3.0;
const double Q_s=-1.0/3.0;


const double GF=1.166389E-5;       //GeV^-2
const double alpha=1.0/137.0359895;
const double PI=3.14159265;

const double proton_mass=0.9383;  //GeV

const double deg_to_rad = 1.0/180.0*PI;


#endif
