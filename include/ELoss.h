#ifndef _ELOSS_H_
#define _ELOSS_H_

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <string.h>

namespace ELOSS
{
    static const double Pi=3.14159265359;
    static const double e=1.60217653e-19;
    static const double Na=6.0221415e23;


    //units, added by Jixie

    /* Common units - Normal : MeV, cm, g, rad */

    const double GeV = 1000.0;
    const double MeV = GeV / 1000.0;
    const double KeV = MeV / 1000.0;
    const double eV  = KeV / 1000.0;

    const double m  = 100.0;
    const double cm = m / 100.0;
    const double mm = m / 1000.0;
    const double um = mm / 1000.0;
    const double nm = um / 1000.0;
    const double pm = nm / 1000.0;
    const double fm = pm / 1000.0;
    const double Km  = 1000.0 * m;
    const double cm2 = cm * cm;
    const double cm3 = cm * cm * cm;

    const double inch = 2.54 *cm;
    const double mil = inch / 1000.0;

    const double barn = 1e-28 * m * m;
    const double milibarn = barn / 1000.0;
    const double microbarn = milibarn / 1000.0;
    const double nanobarn = microbarn / 1000.0;

    const double Kg = 1000.0;
    const double g  = Kg / 1000.0;
    const double mg = g / 1000.0;

    const double rad = 1.0;
    const double mrad = rad / 1000.0;
    const double deg = Pi / 180.0;

    const double percent = 1 / 100.0;


    /* Physical constants - With units */

    const double alpha = 1.0 / 137.035999679; //fine structure static constant
    const double r_e = 2.8179402894e-15*m;
    const double M_e = 0.510998910*MeV;
    const double M_p = 938.272013*MeV;
    const double M_He3 = 2.808391283*GeV;

    const double a2mev=931.49404 * MeV;   //atomic mass unit, transfer A to MeV
    const double GeVSquare2MiliBarn = 0.389379304e6;
    const double GeV2Fm = 0.1973269631;

    struct Material_t
    {
        double Z;
        double A;
        double D;  //density, in g/cm^3 unit
        double L;  //Thickness, in cm unit
        double X0; //radiation length, can be calculated, in g/cm^2 unit
        double M;  //M=A*a2mev;in MeV
        double DL; //density*Thickness,in g/cm^2 unit
        double I;
        std::string name;
        std::string matl_name;//material name
    };

    //////////////////////////////////////////////////////////////////////////
    //set the random seed, if seed=0, will us the time in accuracy of mili second
    void SetRandomSeed(unsigned seed=0);

    //return a random number between 0 and 1
    double Rndm();

    //return a random number in a gaussian distribution
    double RndmGaus();

    ///////////// calculate ionization potential /////////////////////////////
    double GetI(double Z);

    ///////////// calculate radiation length  ////////////////////////////////
    double GetX0(double Z, double A);
    //input: n is number of component
    //     Z[n],A[n],N[n] are the Z|A|number_of_atoms_in_molecular
    //                    lists of components
    double GetX0(int n,double Z[],double A[],int N[]);

    //////////// calculate the density effect: delta /////////////////////////
    //input:
    //E is the energy of incident partical,
    //name and Z are the name and atomic number of the target,
    //para is the parameter which you want to return
    //output: according to para, will return as the following:
    //if(para=="a") return a;
    //else if(para=="m") return m;
    //else if(para=="X0") return X0;
    //else if(para=="X1") return X1;
    //else if(para=="C") return C;
    //else return delta;
    double GetDelta(double E, std::string name,double Z,std::string para);

    void MultipleScattering(double & Theta, double & phi, double Ep, double RadLength);

    //////////////////////////////////////////////////////////////////////////
    //Vivien's MCE97110 C++ energy loss routines
    //return the energy loss in the same unit as E
    double Eloss_IntRad_C(double Ep_init, double Angle, double Mass);
    double Eloss_ExtRad_C(double E, double Z, double Thick_in_Radlen);
    double Eloss_Ion_C(double E, double A, double Z, double Thickness, double Density);

    //Alex's MCE97110 Fortran energy loss routines, tranlated to C by Pengjia
    double Eloss_Ion_F(double E,double A,double Z,double Thick_X_Dens);
    double Eloss_ExtRad_F(double E,double Z,double Thick_in_Radlen);
    double Eloss_IntRad_F(double E0,double theta,double mass);

    //Pengjia's Inoization energy lose routines, the best one we can have
    double Eloss_Ion(double E, Material_t &aMatl);
    double Eloss_Ion_heavy(double E,Material_t &aMatl,double M,int z) ;


    //////////////////////////////////////////////////////////////////////////
    double GetNewE_C(double E, double Z, double A, double L, double D,
        double X0, int WithIntRad=0, double Angle=Pi/2.0);
    double GetNewE_C(double E, Material_t &aMatl,int WithIntRad=0,double ang=Pi/2.0);

    double GetNewE_F(double E, double Z, double A, double L, double D,
        double X0, int WithIntRad=0, double Angle=Pi/2.0);

    double GetNewE_F(double E, Material_t &aMatl,int WithIntRad=0,double Angle=Pi/2.0);

    //Using Pengjia's Ionization energy loss routines
    double GetNewE(double E, Material_t &aMatl,int WithIntRad=0,double Angle=Pi/2.0);

    //Using Pengjia's Ionization energy loss routines for heavy incident particles
    double GetNewE_heavy(double E, Material_t &aMatl,double M,int Z);

    // compares the results among these 3 set of routines
    //input: E0 initial Energy, n is number of materials, aMatl is the pointer to material array
    //       angle is the incident angle, used for internal rad. energy loss
    //output:array E[n], which is the energy of the electron exit each layer of the material
    void CmpELoss(double E0, int nMatl, Material_t *aMatl, int WithIntRad, double angle, double *E);

    ////////////////////////////////////////////////////////////////
    //ToDo: Get the period table of elements
    //Good to have a Z,A table for all elements
    //It will be good to have density table too

    //if I,X0,M not greater than 0, will recalculate them
    void ConstructMaterial(Material_t &aMatl,double Z,double A,double D, double L,char *Matl_name,
        char *name="", double I=0.0,double X0=0.0, double M=0);

    //input: n is number of component
    //     Zlist[n],Alist[n],Nlist[n] are the Z|A|number_of_atoms_in_molecular
    //                    lists of components
    //if I,X0,M not greater than 0, will recalculate them
    void ConstructMaterial(Material_t &aMatl,int n,double Zlist[],double Alist[],int Nlist[],
        double D, double L, char *Matl_name, char *name="", double I=0.0,double X0=0.0, double M=0);

    //if new DL not greater than 0, will get the DL from source material
    //if name is empty, will get the name from source material
    void ConstructMaterial(Material_t &aMatl,const Material_t &srcMatl,double DL=0.0,char *name="");
    ////////////////////////////////////////////////////////////////
  
    //This is the routine you want to use. Let it work as a black box......
    //input: E0 is the initial Energy, n is number of materials, aMatl is the pointer to material array
    //       angle is the incident angle, used for internal rad. energy loss
    //output:array E[n], which is the energy of the electron exit each layer of the material
    
    //calculate Eloss for electron, final energy is in laudao distribution 
    void CalELoss_electron(double E0, int nMatl, Material_t *aMatl, int WithIntRad, double Angle, double *Ef);
    //calculate MEAN Eloss for heavy particles like pion, kaon, proton, alpha ..., , final energy is not in laudao distribution
    void CalELoss_heavy(double E0, int nMatl, Material_t *aMatl, double M, int Z, double *Ef);

}

#endif //_ELOSS_H_
