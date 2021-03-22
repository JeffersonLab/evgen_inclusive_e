//All the constants used in this file are from Leo book or
//Atomic data and nuclear data table 30, 261-271 (1984)
//By R.M. Sternheimer
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "ELoss.h"

using namespace std;

namespace ELOSS
{
    void SetRandomSeed(unsigned seed)
    {
        /*using time now + clock as the seed, whcih is as acurate as mili second*/
        /*If user want to user the same seed every time, just provide a positive number*/
        if(seed==0) seed=unsigned(time(0))+unsigned(clock());
        srand(seed);
    }

    //return the random number between 0 and 1
    double Rndm()
    {
        return double(rand())/double(RAND_MAX);
    }

    double RndmGaus()
    {
        double rndm = Rndm();
        /* Generate a gaussian random number */
        return (rndm>0.0) ? sqrt(-2.0*log(rndm))*cos(2.0*Pi*Rndm()) : 0.0;
    }


    double GetX0_PDG_2720(double A, double Z)
    {
        double Lrad, Lradp, Fz, X0;

        /* Calculation of the factor b */
        if(Z<=2.0)
        {
            Lrad  = (4.79  - 5.31)  * (Z - 1.0) + 5.31;
            Lradp = (5.621 - 6.144) * (Z - 1.0) + 6.144;
        }
        else if(Z<=3.0)
        {
            Lrad  = (4.74  - 4.79)  * (Z - 2.0) + 4.79;
            Lradp = (5.805 - 5.621) * (Z - 2.0) + 5.621;
        }
        else if(Z<=4.0)
        {
            Lrad  = (4.71  - 4.74)  * (Z - 3.0) + 4.74;
            Lradp = (5.924 - 5.805) * (Z - 3.0) + 5.805;
        }
        else
        {
            Lrad  = log(184.15 * pow(Z , -1.0/3.0));
            Lradp = log(1194.0 * pow(Z , -2.0/3.0));
        }

        double t = pow(alpha*Z,2.0);
        Fz = t * ( 1.0/(1+t) + 0.20206 - 0.0369*t + 0.0083*t*t + -0.002*t*t*t );
        X0 = (716.408 * g/cm2) * A / ( Z*Z*(Lrad-Fz) + Z*Lradp );
        cout<<"X0 * A="<<X0<<endl;
        X0 = (716.408 * g/cm2) / ( Z*Z*(Lrad-Fz) + Z*Lradp );
        cout<<"X0 ="<<X0<<endl;
        return X0;
    }

    //return the radiation length in unit of g/cm2
    //From PDG book, eqation 27.22, p273
    double GetX0(double Z,double A)
    {
        return 716.4*A/(Z*(Z+1)*log(287./sqrt(Z))) * g/cm2;
    }

    //return the radiation length in unit of g/cm2
    //input: n is number of component
    //     Z[n],A[n],N[n] are the Z|A|number_of_atoms_in_molecular
    //                    lists of components
    double GetX0(int n,double Z[],double A[],int N[])
    {
        double Xi,Wi;
        double X_this=0.0, W_this=0.0,WDX_sum=0.0;
        for(int i=0;i<n;i++)
        {
            //rad length of each component
            Xi=716.4*A[i]/(Z[i]*(Z[i]+1)*log(287./sqrt((double)Z[i])));
            // weight of each component in 1 mol of this material
            Wi=double(N[i]*A[i]);
            WDX_sum+=Wi/Xi;
            W_this+=Wi;
        }
        X_this=(WDX_sum>0.0)?W_this/WDX_sum:9999999999.0;
        //cout<<"W_this="<<W_this<<"  WDX_sum="<<WDX_sum<<" X_this="<<X_this<<endl;
        return X_this*g/cm2;
    }


    ///////////// calculate ionization potential /////////////////////////////
    double GetI(double Z)
    {
        /* Ionization potentials - source : NIST */
        double I;
        int Z_up,Z_down;
        double Ioni[] = {
            19.2,  41.8,  40.0,  63.7,  76.0,  78.0,  82.0,  95.0,  115.0, 137.0,
            149.0, 156.0, 166.0, 173.0, 173.0, 180.0, 174.0, 188.0, 190.0, 191.0,
            216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0, 311.0, 322.0, 330.0,
            334.0, 350.0, 347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0, 393.0,
            417.0, 424.0, 428.0, 441.0, 449.0, 470.0, 470.0, 469.0, 488.0, 488.0,
            487.0, 485.0, 491.0, 482.0, 488.0, 491.0, 501.0, 523.0, 535.0, 546.0,
            560.0, 574.0, 580.0, 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0,
            694.0, 705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0,
            810.0, 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0,
            878.0, 890.0
        };
        for(int i=0;i<92;i++)
        { Ioni[i] = Ioni[i]*eV; }
        if(Z < 92.0) //from french's code
        {
            Z_down = (int)Z;
            Z_up   = Z_down + 1;
            I = Ioni[Z_down-1] + (Ioni[Z_up-1] - Ioni[Z_down-1]) * (Z - (double)Z_down);
        }
        else
        {
            I = (9.76 * Z + 58.8 * pow(Z,-1.19))*eV;
        }
        return I;
    }


    /***************************************************************************/
    /*                   Internal Bremsstrahlung energy loss                   */
    /***************************************************************************/
    ///Vivien's MCE97110 C++ energy loss routines
    //return the energy loss in the same unit as Ep_init
    double Eloss_IntRad_C(double Ep_init, double Angle, double Mass)
    {
        if(Angle==0.0) return 0.0;

        double Ep_scat, Q2, nu, dE, R;

        /* Scattered energy */
        Ep_scat = Ep_init  / (1.0 + 2.0 * (Ep_init / Mass) * sin(Angle/2.0) * sin(Angle/2.0));

        Q2 = 2.0 * Ep_init * Ep_scat * (1.0 - cos(Angle));

        /* Transfered energy */
        nu = (2.0 * alpha / Pi) * (log(Q2/(M_e*M_e)) - 1.0);

        /* Compute energy loss */
        R = Rndm();
        if(nu>0.0)
            dE = Ep_scat * pow( R*0.999 , (1.0/nu));
        else
            dE = 0.0;

        /* Test to avoid non-reallistic values */
        /* Even if it should not happen        */
        if(dE>Ep_init) dE = Ep_init;
        if(dE<0.0) dE = 0.0;

        return(dE);
    }


    /***************************************************************************/
    /*                   External Bremsstrahlung energy loss                   */
    /***************************************************************************/
    ///Vivien's MCE97110 C++ energy loss routines
    //return the energy loss in the same unit as E
    double Eloss_ExtRad_C(double E, double Z, double Thick_in_Radlen)
    {
        double b, dE, Lrad, Lradp, R;

        /* Calculation of the factor b */
        if(Z<=2.0)
        {
            Lrad  = (4.79  - 5.31)  * (Z - 1.0) + 5.31;
            Lradp = (5.621 - 6.144) * (Z - 1.0) + 6.144;
        }
        else if(Z<=3.0)
        {
            Lrad  = (4.74  - 4.79)  * (Z - 2.0) + 4.79;
            Lradp = (5.805 - 5.621) * (Z - 2.0) + 5.621;
        }
        else if(Z<=4.0)
        {
            Lrad  = (4.71  - 4.74)  * (Z - 3.0) + 4.74;
            Lradp = (5.924 - 5.805) * (Z - 3.0) + 5.805;
        }
        else
        {
            Lrad  = log(184.15 * pow(Z , -1.0/3.0));
            Lradp = log(1194.0 * pow(Z , -2.0/3.0));
        }

        b = (4.0 / 3.0) * (1.0 + (1.0 / 9.0) * (Z + 1.0) / (Lrad * Z + Lradp));

        /* Compute energy loss */
        R = Rndm();
        dE = E * pow(R*0.999,1.0/(b*Thick_in_Radlen));

        /* Test to avoid non-reallistic values */
        /* Even if it should not happen        */
        if(dE>E) dE = E;
        if(dE<0.0) dE = 0.0;

        return(dE);
    }


    /***************************************************************************/
    /*                         Ionization energy loss                          */
    /***************************************************************************/
    ///Vivien's MCE97110 C++ energy loss routines
    //return the energy loss in the same unit as E
    double Eloss_Ion_C(double E, double A, double Z, double Thickness, double Density)
    {
        int Z_down, Z_up;
        double beta2, gamma2, ksi, I, j, delta, u, dE, phil, X, C, hnup;
        double beta, gamma;
        double K = 0.307075*MeV*cm2/g;

        /* parameters to generate a random number following a landau distribution */
        /* Note for later : define those factors outside the function and create  */
        /*                  a landau random number generator function             */
        double lambda[] = {
            -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,
            1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0,  5.5,  6.0,
            6.5,  7.0,  7.5,  8.0,  8.5,  9.0,  9.5, 10.0, 10.5, 11.0,
            11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0,
            16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0,
            21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0, 25.5, 26.0,
            26.5, 27.0, 27.5, 28.0, 28.5, 29.0, 29.5, 30.0
        };
        double intphi[] = {
            0.0000000, 0.0005000, 0.0065000, 0.0290000, 0.0815000, 0.1515000, 0.2390000, 0.3290000, 0.4140000, 0.4875000,
            0.5500000, 0.6025000, 0.6475000, 0.6875001, 0.7195001, 0.7455001, 0.7670001, 0.7850001, 0.8005001, 0.8145001,
            0.8270001, 0.8385001, 0.8482501, 0.8570001, 0.8650001, 0.8722501, 0.8787501, 0.8847501, 0.8902501, 0.8952501,
            0.8998752, 0.9041252, 0.9080002, 0.9115002, 0.9147502, 0.9177502, 0.9205002, 0.9230002, 0.9253752, 0.9276252,
            0.9297502, 0.9317502, 0.9336252, 0.9353752, 0.9370502, 0.9386502, 0.9402002, 0.9417002, 0.9430752, 0.9443253,
            0.9454502, 0.9464502, 0.9474002, 0.9482752, 0.9491102, 0.9498602, 0.9505252, 0.9511502, 0.9517627, 0.9523627,
            0.9529527, 0.9535328, 0.9541028, 0.9546627, 0.9552028, 0.9557328, 0.9562528, 0.9567627
        };

        /* Ionization potentials - source : NIST */
        double Ioni[] = {
            19.2,  41.8,  40.0,  63.7,  76.0,  78.0,  82.0,  95.0,  115.0, 137.0,
            149.0, 156.0, 166.0, 173.0, 173.0, 180.0, 174.0, 188.0, 190.0, 191.0,
            216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0, 311.0, 322.0, 330.0,
            334.0, 350.0, 347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0, 393.0,
            417.0, 424.0, 428.0, 441.0, 449.0, 470.0, 470.0, 469.0, 488.0, 488.0,
            487.0, 485.0, 491.0, 482.0, 488.0, 491.0, 501.0, 523.0, 535.0, 546.0,
            560.0, 574.0, 580.0, 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0,
            694.0, 705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0,
            810.0, 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0,
            878.0, 890.0
        };
        for(int i=0;i<92;i++)
        { Ioni[i] = Ioni[i]*eV; }

        /* Calculation of various parameters */
        //beta2 = 1.0 - (M_e / E) * (M_e / E); //not goog for low E
        beta2 = 1.0 - pow( 1.0/(1.0+E/M_e), 2);
        beta  = sqrt(beta2);
        gamma2 = 1 / (1 - beta2);
        gamma  = sqrt(gamma2);

        ksi = (K / 2.0) * (Z / A) * (Thickness * Density / beta2);

        /* Density correction factor - delta */
        /* Formula are from the PDG booklet and Tsai */
        j = 0.2000;

        if(Z < 92.0)
        {
            Z_down = (int)Z;
            Z_up   = Z_down + 1;

            I = Ioni[Z_down-1] + (Ioni[Z_up-1] - Ioni[Z_down-1]) * (Z - (double)Z_down);
        }
        else
        {
            I = Z * (9.76 + 58.8 * pow(Z,-1.19))*eV;
        }

        hnup = 28.816*eV*sqrt(cm3) * sqrt(Density * Z / A);
        C = -2.0 * log(I/hnup) - 1.0;

        X = log10(beta*gamma);
        delta = 4.6052 * X + C;


        /* Compute energy loss */
        dE = ksi * ( log(2.0 * M_e * beta2 * gamma2 / I)
            + log(ksi / I)
            + j
            - beta2
            - delta
            );

        /* Get random number following a Landau distribution */
        /* Note : do an independant function af this         */
        u = Rndm()*(intphi[67]-intphi[0])+intphi[0];
        phil = 0.0;
        for(int i=1; i<68; i++)
        {
            if(/*u > intphi[i-1] && */u <= intphi[i])
            {
                phil = (lambda[i] - lambda[i-1]) / (intphi[i] - intphi[i-1]);
                phil = phil * (u - intphi[i-1]) + lambda[i-1];
                break;
            }
        }

        /* Compute energy loss with the landau distribution */
        dE = dE + phil * ksi;

        /* Test to avoid non-reallistic values */
        /* Even if it should not happen        */
        if(dE > E) dE = E;
        if(dE < 0.0) dE = 0.0;

        return(dE);
    }

    ////////////////////// by pengjia /////////////////////////////////////////
    //This routine was new written by Pengjia, It is the best ionization energy loss routine
    //
    double Eloss_Ion(double E, Material_t &aMatl)
    {
        double beta2, gamma2, ksi, delta, u, dE, phil, ln_epsi;
        double beta, gamma;
        double K = 0.307075*MeV*cm2/g;
        /* parameters to generate a random number following a landau distribution */
        /* Note for later : define those factors outside the function and create  */
        /*                  a landau random number generator function             */
        double lambda[] = {
            -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,
             1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0,  5.5,  6.0,
             6.5,  7.0,  7.5,  8.0,  8.5,  9.0,  9.5, 10.0, 10.5, 11.0,
            11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0,
            16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0,
            21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0, 25.5, 26.0,
            26.5, 27.0, 27.5, 28.0, 28.5, 29.0, 29.5, 30.0
        };
        double intphi[] = {
            0.0000000, 0.0005000, 0.0065000, 0.0290000, 0.0815000, 0.1515000, 0.2390000, 0.3290000, 0.4140000, 0.4875000,
            0.5500000, 0.6025000, 0.6475000, 0.6875001, 0.7195001, 0.7455001, 0.7670001, 0.7850001, 0.8005001, 0.8145001,
            0.8270001, 0.8385001, 0.8482501, 0.8570001, 0.8650001, 0.8722501, 0.8787501, 0.8847501, 0.8902501, 0.8952501,
            0.8998752, 0.9041252, 0.9080002, 0.9115002, 0.9147502, 0.9177502, 0.9205002, 0.9230002, 0.9253752, 0.9276252,
            0.9297502, 0.9317502, 0.9336252, 0.9353752, 0.9370502, 0.9386502, 0.9402002, 0.9417002, 0.9430752, 0.9443253,
            0.9454502, 0.9464502, 0.9474002, 0.9482752, 0.9491102, 0.9498602, 0.9505252, 0.9511502, 0.9517627, 0.9523627,
            0.9529527, 0.9535328, 0.9541028, 0.9546627, 0.9552028, 0.9557328, 0.9562528, 0.9567627
        };
        /* Calculation of various parameters */
        //beta2 = 1.0 - (M_e / E) * (M_e / E); //not good for low E
        beta2 = 1.0 - pow( 1.0/(1.0+E/M_e), 2);
        beta  = sqrt(beta2);
        gamma2 = 1 / (1 - beta2);
        gamma  = sqrt(gamma2);
        ksi = (K / 2.0) * (aMatl.Z / aMatl.A) * (aMatl.L * aMatl.D/ beta2);
        ln_epsi=log((1-beta2)*pow(aMatl.I,2)/(2*M_e*beta2))+beta2;

        /* Density correction factor - delta */
        /* Formula are from leo's book*/
        delta = GetDelta(E,aMatl.matl_name,aMatl.Z,"delta");
        /* Compute energy loss */
        dE = ksi*(log(ksi)-ln_epsi+0.198-delta);

        /* Get random number following a Landau distribution */
        /* Note : do an independant function af this         */
        u = Rndm()*(intphi[67]-intphi[0])+intphi[0];
        phil = 0.0;
        for(int i=1; i<68; i++)
        {
            if(/*u > intphi[i-1] && */u <= intphi[i])
            {
                phil = (lambda[i] - lambda[i-1]) / (intphi[i] - intphi[i-1]);
                phil = phil * (u - intphi[i-1]) + lambda[i-1];
                break;
            }
        }

        /* Compute energy loss with the landau distribution */
        dE = dE + phil * ksi;

        /* Test to avoid non-reallistic values */
        /* Even if it should not happen        */
        if(dE > E) dE = E;
        if(dE < 0.0) dE = 0.0;

        return(dE);
    }

    //for mean ionization energy loss of heavy partical
    //m and z are the mass and atomic number of the incident partical
    double Eloss_Ion_heavy(double E,Material_t &aMatl,double M,int z)
    {
        double beta2, gamma2, ksi,delta,s,C,Wmax,dEmean;
        double beta, gamma,eta,tmp1,tmp2;
        double K = 0.307075*MeV*cm2/g;
        //from Leo's book 2.27
        /* Calculation of various parameters */
        beta2 = 1.0 - pow( 1.0/(1.0+E/M), 2);
        beta  = sqrt(beta2);
        gamma2 = 1.0 / (1.0 - beta2);
        gamma  = sqrt(gamma2);
        eta=beta*gamma;
        s=M_e/M;
        Wmax=2*M_e*eta*eta/(1+2*s*sqrt(1+(eta*eta))+s*s);
        delta = GetDelta(E,aMatl.matl_name,aMatl.Z,"delta");
        if(delta==9999)
        {
            delta=0;
            C=(0.422377*pow(eta,-2.0) + 0.0304043*pow(eta,-4.0) - 0.00038106*pow(eta,-6.0))*1.0e-6*pow(aMatl.I/eV,2.0) +
              (3.85019*pow(eta,-2.0) - 0.1667989*pow(eta,-4.0) + 0.00157955*pow(eta,-6.0))*1.0e-9*pow(aMatl.I/eV,3.0);
        }
        C=GetDelta(E,aMatl.matl_name,aMatl.Z,"C");
        ksi=(K/2.0)*aMatl.D*aMatl.L*(aMatl.Z/aMatl.A)*pow((z/beta),2);
        tmp1=log(2.*M_e*eta*eta*Wmax/pow(aMatl.I,2));
        tmp2=tmp1-2.*beta2-delta-2.*C/aMatl.Z;
        dEmean=ksi*tmp2;

        if(dEmean>E) dEmean=E;
        if(dEmean<0.0) dEmean=0.0;

        return dEmean;
    }


    //////////////////// calculate the density effect: delta /////////////////////////////////

    //input:
    //E is the energy of incident partical,
    //name and Z are the name and atomic number of the target,
    //para is the parameter which want to return
    //output: according to para, will return as the following:
    //if(para=="a") return a;
    //else if(para=="m") return m;
    //else if(para=="X0") return X0;
    //else if(para=="X1") return X1;
    //else if(para=="C") return C;
    //else return delta;
    //this routine can run faster than the general GetDelta() since it does not look thru all 98 materials
    double GetDeltaG2P(double E,std::string name,double Z,std::string para)
    {
        //data from Atomic data and nuclear data tables 30,261-271(1984)
        double X,gamma2,be2,eta,delta,X0,X1,C,a,m,delta0;
        int i;
        static const int NComp=3;//how many componants have data in here
        bool findmaterial=false;//judge whether the material has been identified
        bool isconductor=false; //judge whether the material is a conductor

        string material_name2[]={"NH3","BeO","Air"}; //componant


        double X0s2[] = {1.6822,0.0241,1.7418};
        double X1s2[] = {4.1158,2.5846,4.2759};
        double Cs2[] = {9.8763,2.9801,10.5961};
        double as2[] = {0.08315,0.10755,0.10914};
        double ms2[] = {3.6464,3.4927,3.3994};


        for(i=0;i<NComp;i++)
        {//for componant
            if(name==material_name2[i])
            {
                X0=X0s2[i];
                X1=X1s2[i];
                C=-Cs2[i];
                a=as2[i];
                m=ms2[i];
                delta0=0;
                findmaterial=true;
                break;
            }
        }

        if(!findmaterial)
        {
            std::cout<<"Can not identify the material, please add material parameters into this code"<<std::endl;
            return 9999;
        }

        if(para=="a") return a;
        else if(para=="m") return m;
        else if(para=="X0") return X0;
        else if(para=="X1") return X1;
        else if(para=="C") return C;
        else
        {
            be2=1.0-pow( 1.0/(1.0+E/M_e), 2.0);
            gamma2=1.0/(1.0-be2);
            X=log10(sqrt(be2*gamma2));
            eta=sqrt(be2*gamma2);


            if(X>X1) delta=4.6052*X+C;
            else if(X<X1&&X>X0)  delta=4.6052*X+a*pow((X1-X),m)+C;
            else delta=0;

            //    std::cout<<name<<"  "<<X0<<"  "<<X1<<"  "<<C<<"  "<<a<<"  "<<m<<"  "<<std::endl;

            return delta;
        }
    }


    //////////////////// calculate the density effect: delta /////////////////////////////////
    //input:
    //E is the energy of incident partical,
    //name and Z are the name and atomic number of the target,
    //para is the parameter which want to return
    //output: according to para, will return as the following:
    //if(para=="a") return a;
    //else if(para=="m") return m;
    //else if(para=="X0") return X0;
    //else if(para=="X1") return X1;
    //else if(para=="C") return C;
    //else return delta;

    double GetDelta(double E,std::string name,double Z,std::string para)
    {
        //data from Atomic data and nuclear data tables 30,261-271(1984)
        double X,gamma2,be2,eta,delta,X0,X1,C,a,m,delta0;
        int i;
        static const int NComp=3;//how many componants have data in here
        bool findmaterial=false;//judge whether the material has been identified
        bool isconductor=false; //judge whether the material is a conductor

        std::string material_name[98] = {
            "H","LH","He","Li","Be","B","C_2.3","C_2.0","C_1.7","N",
            "O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
            "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co",
            "Ni","Cu","Zn","Ca","Ge","As","Se","Br","Kr","Rb",
            "Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
            "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La",
            "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho",
            "Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir",
            "Pt","Au","Hg","Tl","Pb","Bi","Po","Rn","Ra","Ac",
            "Th","Pa","U","Np","Pu","Am","Cm","Bk"};
            //different density of C has different parameters,also for liquid H

        string material_name2[]={"NH3","BeO","Air"}; //componant

        double Zs[98] = {
            1,1,2,3,4,5,6,6,6,7,8,9,
            10,11,12,13,14,15,16,17,18,19,
            20,21,22,23,24,25,26,27,28,29,
            30,31,32,33,34,35,36,37,38,39,
            40,41,42,43,44,45,46,47,48,49,
            50,51,52,53,54,55,56,57,58,59,
            60,61,62,63,64,65,66,67,68,69,
            70,71,72,73,74,75,76,77,78,79,
            80,81,82,83,84,86,88,89,90,91,
            92,93,94,95,96,97};

        double Z_conductor[76] = {
            3,4,6,11,12,13,19,20,21,22,
            23,24,25,26,27,28,29,30,31,32,
            37,38,39,40,41,42,43,44,45,46,
            47,48,49,50,51,55,56,57,58,59,
            60,61,62,63,64,65,66,67,68,69,
            70,71,72,73,74,75,76,77,78,79,
            80,81,82,83,84,87,88,89,90,91,
            92,93,94,95,96,97};//conductors

        double X0s[98] = {
            1.8639,0.4759,2.2017,0.1304,0.0592,0.0305,-0.0178,-0.0351, 0.048,1.7378,
            1.7541,1.8433,2.0735, 0.288,0.1499,0.1708, 0.2014, 0.1696, 0.158,1.5555,
            1.7635,0.3851,0.3228, 0.164,0.0957,0.0691,0.034,0.0447,-0.0012,-0.0187,
            -0.0566,-0.0254,0.0049,0.2267,0.3376,0.1767,0.2258,1.5262,1.7158,0.5737,
            0.4585,0.3608,0.2957,0.1785,0.2267,0.0949,0.0599,0.0576,0.0563,0.0657,
            0.1281,0.2406,0.2879,0.3189,0.3296,0.0549,1.563,0.5473,0.419,0.3161,
            0.2713,0.2333,0.1984,0.1627,0.152,0.1888,0.1058,0.0947,0.0822,0.0761,
            0.0648,0.0812,0.1199,0.156,0.1965,0.2117,0.2167,0.0559,0.0891,0.0819,
            0.1484,0.2021,0.2756,0.3491,0.3776,0.4152,0.4267,1.5368,0.5991,0.4559,
            0.4202,0.3144,0.226,0.1869,0.1557,0.2274,0.2484,0.2378};
        double X0s2[]  = {1.6822,0.0241,1.7418};

        double X1s[98] = {
            3.2718,1.9215,3.6122,1.6397,1.6922,1.9688,2.3415,2.486,2.5387,4.1323,
            4.3213,4.4096,4.6421,3.1962,3.0668,3.0127,2.8715,2.7815,2.7159,4.2996,
            4.4855,3.1724,3.1191,3.0593,3.0386,3.0322,3.0451,3.1074,3.1531,3.179,
            3.1851,3.2792,3.3668,3.5434,3.6096,3.5702,3.6264,4.9899,5.0748,3.7995,
            3.6778,3.5542,3.489,3.2201,3.2784,3.1253,3.0834,3.1069,3.0555,3.1074,
            3.1667,3.2032,3.2959,3.3489,3.4418,3.2596,4.7371,3.5914,3.4547,3.3293,
            3.3432,3.2773,3.3063,3.3199,3.346,3.4633,3.3932,3.4224,3.4474,3.4782,
            3.4922,3.5085,3.6246,3.5218,3.4337,3.4805,3.496,3.4845,3.5414,3.548,
            3.6212,3.6979,3.7275,3.8044,3.8073,3.8248,3.8293,4.9889,3.9428,3.7966,
            3.7681,3.5079,3.3721,3.369,3.3981,3.5021,3.516,3.5186};
        double X1s2[]  = {4.1158,2.5846,4.2759};

        double Cs[98]  = {
            9.5835,3.2632,11.1393,3.1221,2.7847,2.8477,2.868,2.9925,3.155,10.54,
            10.7004,10.9653,11.9041,5.0526,4.5297,4.2395,4.4351,4.5214,4.6659,
            11.1421,11.948,5.6423,5.0396,4.6949,4.445,4.2659,4.1781,4.2702,4.2911,
            4.2601,4.3115,4.419,4.6906,4.9353,5.1411,5.051,5.321,11.7307,12.5115,
            6.4776,5.9867,5.4801,5.1774,5.0141,4.8793,4.7769,4.7694,4.8008,4.9358,
            5.063,5.2727,5.5211,5.534,5.6241,5.7131,5.9488,12.7281,6.9135,6.3153,
            5.785,5.7837,5.8096,5.829,5.8224,5.8597,6.2278,5.8738,5.9045,5.9183,
            5.9587,5.9521,5.9677,6.3325,5.9785,5.7139,5.5262,5.4059,5.3445,5.3083,
            5.3418,5.4732,5.5747,5.9605,6.1365,6.2018,6.3505,6.4003,13.2839,7.0452,
            6.3742,6.2473,6.0327,5.8694,5.8149,5.8748,6.2813,6.3097,6.2912};
        double Cs2[]   ={9.8763,2.9801,10.5961};

        double as[98]  = {
            0.14092,0.13483,0.13443,0.95136,0.80392,0.56224,0.26142,0.2024,0.20762,0.15349,
            0.11778,0.11083,0.08064,0.07772,0.08163,0.08024,0.14921,0.2361,0.33992,0.19849,
            0.19714,0.19827,0.15643,0.15754,0.15662,0.15436,0.15419,0.14973,0.1468,0.14474,
            0.16496,0.14339,0.14714,0.0944,0.07188,0.06633,0.06568,0.06335,0.07446,0.07261,
            0.07165,0.07138,0.07177,0.13863,0.10525,0.16572,0.19342,0.19205,0.24178,0.24585,
            0.24609,0.23879,0.18689,0.16652,0.13815,0.23766,0.23314,0.18233,0.18268,0.18591,
            0.18885,0.23265,0.2353,0.2428,0.24698,0.24448,0.25109,0.24453,0.24665,0.24638,
            0.24823,0.24889,0.25295,0.24033,0.22918,0.17798,0.15509,0.15184,0.12751,0.1269,
            0.11128,0.09756,0.11014,0.09455,0.09359,0.0941,0.09282,0.20798,0.08804,0.08567,
            0.08655,0.1477,0.19677,0.19741,0.20419,0.20308,0.20257,0.20192};
        double as2[]={0.08315,0.10755,0.10914};

        double ms[98]  = {
            5.7273,5.6249,5.8347,2.4993,2.4339,2.4512,2.8697,3.0036,2.9532,3.2125,
            3.2913,3.2962,3.5771,3.6452,3.6166,3.6345,3.2546,2.9158,2.6456,2.9702,
            2.9618,2.9233,3.0745,3.0517,3.0302,3.0163,2.9896,2.9796,2.9632,2.9502,
            2.843,2.9044,2.8652,3.1314,3.3306,3.4176,3.4317,3.467,3.4051,3.4177,
            3.4435,3.4585,3.4533,3.093,3.2549,2.9738,2.8707,2.8633,2.7239,2.6899,
            2.6772,2.7144,2.8576,2.9319,3.0354,2.7276,2.7414,2.8866,2.8906,2.8828,
            2.8592,2.7331,2.705,2.6674,2.6403,2.6245,2.5977,2.6056,2.5849,2.5726,
            2.5573,2.5469,2.5141,2.5643,2.6155,2.7623,2.8447,2.8627,2.9608,2.9658,
            3.0417,3.1101,3.0519,3.145,3.1608,3.1671,3.183,2.7409,3.2454,3.2683,
            3.261,2.9845,2.8171,2.8082,2.7679,2.7615,2.7579,2.756};
        double ms2[]   ={3.6464,3.4927,3.3994};

        double delta0s[98]={
            0.00,0.00,0.00,0.14,0.14,0.14,0.12,0.10,0.14,0.00,
            0.00,0.00,0.00,0.08,0.08,0.12,0.14,0.14,0.14,0.00,
            0.00,0.10,0.14,0.10,0.12,0.14,0.14,0.14,0.12,0.12,
            0.10,0.08,0.08,0.14,0.14,0.08,0.10,0.00,0.00,0.14,
            0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,
            0.14,0.14,0.14,0.14,0.14,0.00,0.00,0.14,0.14,0.14,
            0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,
            0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.08,0.10,0.10,
            0.12,0.14,0.14,0.14,0.14,0.14,0.14,0.00,0.14,0.14,
            0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14};

        for(i=0;i<98;i++)
        {
            if(name==material_name[i])
            {
                X0=X0s[i];
                X1=X1s[i];
                C=-Cs[i];
                a=as[i];
                m=ms[i];
                delta0=delta0s[i];
                findmaterial=true;
                break;
            }
        }

        for(i=0;i<NComp;i++)
        {//for componant
            if(name==material_name2[i])
            {
                X0=X0s2[i];
                X1=X1s2[i];
                C=-Cs2[i];
                a=as2[i];
                m=ms2[i];
                delta0=0;
                findmaterial=true;
                break;
            }
        }

        if(!findmaterial)
        { //search by Z
            for(i=0;i<98;i++)
            {
                if(Z==Zs[i])
                {
                    X0=X0s[i];
                    X1=X1s[i];
                    C=-Cs[i];
                    a=as[i];
                    m=ms[i];
                    delta0=delta0s[i];
                    findmaterial=true;
                    break;
                }
            }
        }

        if(!findmaterial)
        {
            std::cout<<"Can not identify the material, please add material parameters into this code"<<std::endl;
            return 9999;
        }

        if(para=="a") return a;
        else if(para=="m") return m;
        else if(para=="X0") return X0;
        else if(para=="X1") return X1;
        else if(para=="C") return C;
        else
        {
            be2=1.0-pow( 1.0/(1.0+E/M_e), 2.0);
            gamma2=1.0/(1.0-be2);
            X=log10(sqrt(be2*gamma2));
            eta=sqrt(be2*gamma2);
            for(i=0;i<76;i++)
            {
                if(Z==Z_conductor[i]) { isconductor=true; break; }
            }
            if(X>X1)
            {
                delta=4.6052*X+C;
            }
            else if(X<X1&&X>X0)
            {
                delta=4.6052*X+a*pow((X1-X),m)+C;
            }
            else
            {
                if(isconductor) delta=pow(10,(2*(X-X0)))*delta0;
                else delta=0;
            }
            //    std::cout<<name<<"  "<<X0<<"  "<<X1<<"  "<<C<<"  "<<a<<"  "<<m<<"  "<<std::endl;

            return delta;
        }
    }

    /*****************************************************************************/
    /*                            Multiple Scattering                            */
    /*****************************************************************************/

    void MultipleScattering(double &theta, double &phi, double Ep, double Thick_in_Radlen)
    {
        double tmp, angle_sigma, angle_scat, phi_scat, theta_scat;
        double Es = 13.6*MeV;

        /* Formula are from the PDG booklet     */
        /* Average width of the angle dipersion */
        if(Thick_in_Radlen>0.0)
            angle_sigma = Es * sqrt(Thick_in_Radlen) / Ep * (1.0 + 0.038 * log(Thick_in_Radlen));
        else
            angle_sigma = 0.0;

        tmp = RndmGaus();
        angle_scat = tmp * angle_sigma;

        tmp = Rndm();
        phi_scat   = angle_scat * sin(2.0 * Pi * tmp);
        theta_scat = angle_scat * cos(2.0 * Pi * tmp);

        /* New theta and phi */
        theta = theta + theta_scat;
        phi   = phi + phi_scat;
    }

    /////////////////////////////////////////////////////////////////////////
    //The following routine is translated from Alex's MCE97110 fortran code by Pengjia
    double Eloss_IntRad_F(double E0,double theta,double mass)
    {
        if(E0<=0.0) return 0.0;
        if(theta==0.0) return 0.0;

        double EE,dEE,qsq,nu,tmp;
        EE=E0/(1.+2.*E0/mass*pow(sin(theta/2.),2));
        qsq=2.*EE*E0*(1.-cos(theta));
        nu=alpha/Pi*(log(qsq/M_e/M_e)-1.);
        tmp=Rndm()*0.999;
        dEE=EE*pow(tmp,(1./nu));

        if(dEE>E0) dEE=E0;
        return dEE;
    }

    //The following routine is translated from Alex's MCE97110 fortran code by Pengjia
    double Eloss_ExtRad_F(double E,double z,double Thick_in_Radlen)
    {
        if(E<=0.0) return 0.0;

        double b,tmp;
        double dE;

        b=4./3.*(1.+(z+1.)/(9.*(log(1194.*pow(z,(-2./3.)))+z*log(184.15*pow(z,(-1./3.))))));
        tmp=Rndm()*0.999;
        dE=E*pow(tmp,(1/(b*Thick_in_Radlen)));

        if(dE>E) dE=E;

        return dE;
    }

    //The following routine is translated from Alex's MCE97110 fortran code by Pengjia
    double Eloss_Ion_F(double E,double A,double Z,double Thick_X_Dens)
    {
        if (E<=0.0) return 0.0;

        int i;
        double tmp;
        double be2,be2mo,dE,dEprob,phil;
        double ksi,lneps;
        double ld[2][68]={
            {-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-.5,.0,.5,1.0,1.5,2.,2.5,3.0,3.5,4.0,
            4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0 ,11.5,12.0,
            12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0,16.5,17.,17.5,18.,18.5,19.,
            19.5,20.,20.5,21.,21.5,22.,22.5,23.,23.5,24.0,24.5,25.0,25.5,26.,
            26.5,27.,27.5,28.,28.5,29.,29.5,30.},
            {0.0,5.0e-04,6.5e-03,2.90e-02,8.15e-02,.1515,.239,.329,.414,
            .4875,.55,.6025,.6475,.6875001,.7195001,.7455001,.7670001,.7850001,.8005001,
            .8145001,.8270001,.8385001,.8482501,.8570001,.8650001,.8722501,.8787501,
            .8847501,.8902501,.8952501,.8998752,.9041252,.9080002,.9115002,.9147502,
            .9177502,.9205002,.9230002,.9253752,.9276252,.9297502,.9317502,.9336252,
            .9353752,.9370502,.9386502,.9402002,.9417002,.9430752,.9443253,.9454502,
            .9464502,.9474002,.9482752,.9491102,.9498602,.9505252,.9511502,.9517627,
            .9523627,.9529527,.9535328,.9541028,.9546627,.9552028,.9557328,.9562528,.9567627}
        };

        //    be2=1.-pow((M_e/E),2);
        //    be2mo=pow((M_e/E),2);
        be2=1.0-pow(1.0/(1.0+E/M_e),2);
        //be2=be2Ca(E);
        be2mo=1-be2;
        ksi=0.154/be2*Z/A*Thick_X_Dens;
        lneps=log(be2mo*pow(0.0000135*Z,2)/2./M_e/be2)+be2;
        dEprob=ksi*(log(ksi)-lneps+0.37);
        tmp=Rndm()*(ld[1][67]-ld[1][0])+ld[1][0];
        for(i=1;i<68;i++)
        {
            if(/*tmp>ld[1][i-1] && */tmp<=ld[1][i])
            {
                phil=(ld[0][i]-ld[0][i-1])/(ld[1][i]-ld[1][i-1])*(tmp-ld[1][i-1])+ld[0][i-1];
                break;
            }
        }

        dE=dEprob+0.06+phil*ksi;

        //if(dE<=-0.||dE>=E)
        //{
        //    cout<<"Get A wrong number in Eloss_Ion_F"<<endl;
        //    cout<<"E="<<E<<"  dE="<<dE<<"  dEprob="<<dEprob<<"  phil="<<phil
        //        <<"  ksi="<<ksi<<endl;
        //}
        if(dE>E) dE=E;
        if(dE<0) dE=0.0;

        return dE;
    }
    //The above routines are translated from the MCE97110 fortran code
    ///////////////////////////////////////////////////////////////////////////


    //The routine use Vivien's MCE97110 C++ code
    double GetNewE_C(double E, double Z, double A, double L, double D,
        double X0, int WithIntRad, double Angle)
    {
        double dEextr=Eloss_ExtRad_C(E,Z,L*D/X0);
        //cout<<"dEextr="<<dEextr<<endl;
        if(E>dEextr && dEextr>=0) E-=dEextr;
        else return 0.0;

        if(WithIntRad)
        {
            double amu=931.49404; //atomic mass unit
            double mass=A*amu;
            double dEintr=Eloss_IntRad_C(E,Angle,mass);
            //cout<<"dEintr="<<dEintr<<endl;
            if(E>dEintr && dEintr>=0) E-=dEintr;
            else return 0.0;
        }

        double dEi=Eloss_Ion_C(E,A,Z,L,D);
        //cout<<"dEi="<<dEi<<endl;
        E=(E>dEi && dEi>=0)?E-dEi:0.0;
        return E;
    }

    //The routine use Viven's MCE97110 C++ code
    double GetNewE_C(double E, Material_t &aMatl,int WithIntRad, double Angle)
    {
        double dEextr=Eloss_ExtRad_C(E,aMatl.Z,aMatl.DL/aMatl.X0);
        //cout<<"dEextr="<<dEextr<<endl;
        if(E>dEextr && dEextr>=0) E-=dEextr;
        else return 0.0;
        if(WithIntRad)
        {
            double dEintr=Eloss_IntRad_C(E,Angle,aMatl.M);
            //cout<<"dEintr="<<dEintr<<endl;
            if(E>dEintr && dEintr>=0) E-=dEintr;
            else return 0.0;
        }

        double dEi=Eloss_Ion_C(E,aMatl.A,aMatl.Z,aMatl.L,aMatl.D);
        //cout<<"dEi="<<dEi<<endl;
        E=(E>dEi && dEi>=0)?E-dEi:0.0;
        return E;
    }

    //The routine use the translated Alex's MCE97110 Fortran code
    double GetNewE_F(double E, double Z, double A,  double L, double D,
        double X0,int WithIntRad, double Angle)
    {
        double dEextr=Eloss_ExtRad_F(E,Z,L*D/X0);
        if(E>dEextr && dEextr>=0) E-=dEextr;
        else
        {
            cout<<"E="<<E<<"\t dEextr="<<dEextr<<endl;
            return 0.0;
        }

        if(WithIntRad)
        {
            double amu=931.49404; //atomic mass unit
            double mass=A*amu;
            double dEintr=Eloss_IntRad_F(E,Angle,mass);
            //cout<<"dEintr="<<dEintr<<endl;
            if(E>dEintr && dEintr>=0) E-=dEintr;
            else return 0.0;
        }

        double dEi=Eloss_Ion_F(E,A,Z,L*D);
        if(E>=dEi && dEi>=0) E-=dEi;
        else
        {
            cout<<"E="<<E<<"\t dEi="<<dEi<<endl;
            return 0.0;
        }
        return E;
    }

    //The routine use the translated Alex's MCE97110 Fortran code
    double GetNewE_F(double E, Material_t &aMatl,int WithIntRad,double Angle)
    {
        double dEextr=Eloss_ExtRad_F(E,aMatl.Z,aMatl.DL/aMatl.X0);
        //cout<<"dEextr="<<dEextr<<endl;
        if(E>dEextr && dEextr>=0) E-=dEextr;
        else return 0.0;

        if(WithIntRad)
        {
            double dEintr=Eloss_IntRad_F(E,Angle,aMatl.M);
            //cout<<"dEintr="<<dEintr<<endl;
            if(E>dEintr && dEintr>=0) E-=dEintr;
            else return 0.0;
        }

        double dEi=Eloss_Ion_F(E,aMatl.A,aMatl.Z,aMatl.DL);
        //cout<<"dEi="<<dEi<<endl;
        E=(E>dEi && dEi>=0)?E-dEi:0.0;
        return E;
    }

    ////////////////////// pengjia ///////////////////////////////////////////////
    //This routine used Pengjia's code
    double GetNewE(double E, Material_t &aMatl,int WithIntRad,double Angle)
    {
        //double dEextr=Eloss_ExtRad_F(E,aMatl.Z,aMatl.DL/aMatl.X0);
        double dEextr=Eloss_ExtRad_C(E,aMatl.Z,aMatl.DL/aMatl.X0);
        //cout<<"dEextr="<<dEextr<<endl;
        if(E>dEextr && dEextr>=0) E-=dEextr;
        else return 0.0;
	
        if(WithIntRad)
        {
            double dEintr=Eloss_IntRad_F(E,Angle,aMatl.M);
            //cout<<"dEintr="<<dEintr<<endl;
            if(E>dEintr && dEintr>=0) E-=dEintr;
            else return 0.0;
        }

        double dEi=Eloss_Ion(E,aMatl);
        //double dEi=Eloss_Ion_F(E,aMatl.A,aMatl.Z,aMatl.DL);
        //cout<<"dEi="<<dEi<<endl;
        E=(E>dEi && dEi>=0)?E-dEi:0.0;
        return E;
    }

    //Get the energy after ionization for heay particles like proton, Helium and Kaons
    //M and Z are the mass and atomic number of the incident partical
    double GetNewE_heavy(double E, Material_t &aMatl,double M,int Z)
    {
        double dEi=Eloss_Ion_heavy(E,aMatl,M,Z);
        E=(E>dEi && dEi>=0)?E-dEi:0.0;
        return E;
    }


    //this is used for building single atom molecula
    void ConstructMaterial(Material_t &aMatl,double Z,double A,double D, double L,char *matl_name,
        char *name,double I,double X0, double M)
    {
        aMatl.Z=Z;
        aMatl.A=A;
        aMatl.D=D;
        aMatl.L=L;
        aMatl.DL=D*L;
        aMatl.I=(I<=0.0&&Z>0)?GetI(Z):I;
        aMatl.M=(M<=0.0) ? a2mev*A : M;
        aMatl.X0=(X0<=0.0 && Z>0) ? GetX0(Z,A) : X0;
        aMatl.name=name;
        aMatl.matl_name=matl_name;
    }


    //this is used for building multi-atoms molecula
    void ConstructMaterial(Material_t &aMatl,int n,double Zlist[],double Alist[],int Nlist[],
        double D, double L,char *matl_name,char *name,double I,double X0, double M)
    {
        //calculate average Z and average A
        double Wi, W_this=0.0,WXA_sum=0.0, WXZ_sum=0.0;
        for(int i=0;i<n;i++)
        {
            // weight of each component in 1 mol of this material
            Wi=double(Nlist[i])*Alist[i];
            WXA_sum+=Wi*Alist[i];
            WXZ_sum+=Wi*(double)Zlist[i];
            W_this+=Wi;
        }
        if(W_this<=0)
        {
            std::cout<<"Error in ConstructMaterial(), invalid input arguments"<<std::endl;
            return;
        }
        double Z_av=WXZ_sum/W_this;
        double A_av=WXA_sum/W_this;

        aMatl.Z=Z_av;
        aMatl.A=A_av;
        aMatl.D=D;
        aMatl.L=L;
        aMatl.I=(I<=0.0&&Z_av>0)?GetI(Z_av):I;
        aMatl.DL=D*L;
        aMatl.M=(M<=0.0) ? a2mev*W_this : M;
        aMatl.X0=(X0<=0.0 && Z_av>0) ? GetX0(n,Zlist,Alist,Nlist) : X0;
        aMatl.name=name;
        aMatl.matl_name=matl_name;
    }

    //build a new material base on the given one, reset the thickness if given L>0
    void ConstructMaterial(Material_t &aMatl,const Material_t &srcMatl,double L, char *name)
    {
        aMatl.Z=srcMatl.Z;
        aMatl.A=srcMatl.A;
        aMatl.D=srcMatl.D;
        aMatl.L=(L>0)?L:srcMatl.L;
        aMatl.X0=srcMatl.X0;
        aMatl.M=srcMatl.M;
        aMatl.DL=aMatl.D*aMatl.L;
        aMatl.I=srcMatl.I;
        aMatl.name=(strlen(name)>0)?name:srcMatl.name;
        aMatl.matl_name=srcMatl.matl_name;
    }

    //This is an example of how to use this code. It also compares the results among these
    //3 set of routines
    //input: E0 is the initial Energy, n is number of materials, aMatl is the pointer to material array
    //       angle is the incident angle, used for internal rad. energy loss
    //output:array E[n], which is the energy of the electron exit each layer of the material
    void CmpELoss(double E0, int nMatl, Material_t *aMatl, int WithIntRad, double Angle, double *Ee)
    {
        double *tmpEe_C,*tmpEe_F;
        tmpEe_C=new double [nMatl];
        tmpEe_F=new double [nMatl];
        double diff_C,diff_F;
        for(int i=0;i<nMatl;i++)
        {
            Ee[i]=(i>0)?Ee[i-1]:E0;
            if(aMatl[i].A > 0 && Ee[i]>0.0)
            {
                //tmpEe_C[i]=ELOSS::GetNewE_C(E[i],aMatl[i].Z,aMatl[i].A,aMatl[i].D,aMatl[i].L,aMatl[i].X0,WithIntRad,Angle);
                //tmpEe_F[i]=ELOSS::GetNewE_F(E[i],aMatl[i].Z,aMatl[i].A,aMatl[i].D,aMatl[i].L,aMatl[i].X0,WithIntRad,Angle);
                tmpEe_C[i]=ELOSS::GetNewE_C(Ee[i],aMatl[i],WithIntRad,Angle);
                tmpEe_F[i]=ELOSS::GetNewE_F(Ee[i],aMatl[i],WithIntRad,Angle);
                //Ee[i]=GetNewE_C(E[i],aMatl[i],WithIntRad,Angle);
                Ee[i]=GetNewE(Ee[i],aMatl[i],WithIntRad,Angle);
		diff_C=Ee[i]-tmpEe_C[i];
                diff_F=Ee[i]-tmpEe_F[i];
                cout<<"diff_C="<<diff_C<<"  diff_F="<<diff_F<<"  Ee[i]="<<Ee[i]<<endl;
            }
        }
        delete [] tmpEe_C;
        delete [] tmpEe_F;
    }
   
    ////////////////////////////////////////////////////////////////
    //This is the routine you want to use. Let it work as a black box......
    //input: E0 is the initial Energy, n is number of materials, aMatl is the pointer to material array
    //       angle is the incident angle, used for internal rad. energy loss
    //output:array E[n], which is the energy of the electron exit each layer of the material
    
    //calculate Eloss for electron, final energy is in laudao distribution 
    void CalELoss_electron(double E0, int nMatl, Material_t *aMatl, int WithIntRad, double Angle, double *Ef)
    {
        for(int i=0;i<nMatl;i++)
        {
            Ef[i]=(i>0)?Ef[i-1]:E0;
            if(aMatl[i].A > 0 && Ef[i]>0.0)
            {
                Ef[i]=GetNewE(Ef[i],aMatl[i],WithIntRad,Angle);
            }
        }
    } 
    
    //calculate MEAN Eloss for heavy particles like pion, kaon, proton, alpha ..., , final energy is not in laudao distribution
    void CalELoss_heavy(double E0, int nMatl, Material_t *aMatl, double M, int Z, double *Ef)
    {
        for(int i=0;i<nMatl;i++)
        {
            Ef[i]=(i>0)?Ef[i-1]:E0;
            if(aMatl[i].A > 0 && Ef[i]>0.0)
            {
                Ef[i]=GetNewE_heavy(Ef[i],aMatl[i],M,Z);
            }
        }
    }

}  //end of namespace
