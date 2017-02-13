
#include "math.h"
#include <cmath>
#include "iostream"
#include "christy_bosted_inelastic_QE.h"

using namespace std;


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

/*

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

*/



// the real code to get F1 and F2 calculation...
int F1F2IN09(int Z, int IA, double qsq, double wsq, double &F1, double &F2, double &Rc){  //Rc removed from the Frotran code
    // q2 and W are using GeV 
    
    /*--------------------------------------------------------------------
     Fit to inelastic cross sections for A(e,e')X
     valid for all W<3 GeV and all Q2<10 GeV2

     Inputs: Z, A (real*8) are Z and A of nucleus
             (use Z=0., A=1. to get free neutron)
             Qsq (real*8) is 4-vector momentum transfer squared (positive in
                         chosen metric)
             Wsq (real*8) is invarinat mass squared of final state calculated
                         assuming electron scattered from a free proton

     outputs: F1, F2 (real*8) are structure functions per nucleus
     Version of 10/20/2006 P. Bosted
    --------------------------------------------------------------------*/

    double W, x;
    double nu,siglp,sigtp,F1pp,F1dp;
    double W1,W2,sigt,sigl,F1d, F1p,qv;
    double DW2DPF,Wsqp,pf,kf,Es,dw2des,Fyuse;
    double x4, emcfac;
    int ism;
    double PM = 0.938;

    // This is for exp(-xx**2/2.), from teste.f
    double XXp[99] = {
        -3.000,-2.939,-2.878,-2.816,-2.755,-2.694,-2.633,-2.571,-2.510,
        -2.449,-2.388,-2.327,-2.265,-2.204,-2.143,-2.082,-2.020,-1.959,
        -1.898,-1.837,-1.776,-1.714,-1.653,-1.592,-1.531,-1.469,-1.408,
        -1.347,-1.286,-1.224,-1.163,-1.102,-1.041,-0.980,-0.918,-0.857,
        -0.796,-0.735,-0.673,-0.612,-0.551,-0.490,-0.429,-0.367,-0.306,
        -0.245,-0.184,-0.122,-0.061, 0.000, 0.061, 0.122, 0.184, 0.245,
        0.306, 0.367, 0.429, 0.490, 0.551, 0.612, 0.673, 0.735, 0.796,
        0.857, 0.918, 0.980, 1.041, 1.102, 1.163, 1.224, 1.286, 1.347,
        1.408, 1.469, 1.531, 1.592, 1.653, 1.714, 1.776, 1.837, 1.898,
        1.959, 2.020, 2.082, 2.143, 2.204, 2.265, 2.327, 2.388, 2.449,
        2.510, 2.571, 2.633, 2.694, 2.755, 2.816, 2.878, 2.939, 3.000
    };
    // these are 100x bigger for convenience
    double fyp[99] = {
        0.0272,0.0326,0.0390,0.0464,0.0551,0.0651,0.0766,0.0898,0.1049,
        0.1221,0.1416,0.1636,0.1883,0.2159,0.2466,0.2807,0.3182,0.3595,
        0.4045,0.4535,0.5066,0.5637,0.6249,0.6901,0.7593,0.8324,0.9090,
        0.9890,1.0720,1.1577,1.2454,1.3349,1.4254,1.5163,1.6070,1.6968,
        1.7849,1.8705,1.9529,2.0313,2.1049,2.1731,2.2350,2.2901,2.3379,
        2.3776,2.4090,2.4317,2.4454,2.4500,2.4454,2.4317,2.4090,2.3776,
        2.3379,2.2901,2.2350,2.1731,2.1049,2.0313,1.9529,1.8705,1.7849,
        1.6968,1.6070,1.5163,1.4254,1.3349,1.2454,1.1577,1.0720,0.9890,
        0.9090,0.8324,0.7593,0.6901,0.6249,0.5637,0.5066,0.4535,0.4045,
        0.3595,0.3182,0.2807,0.2466,0.2159,0.1883,0.1636,0.1416,0.1221,
        0.1049,0.0898,0.0766,0.0651,0.0551,0.0464,0.0390,0.0326,0.0272
    };

    // deuteron fit parameters
    double xvald0[50] = {
        0.1964E+01, 0.1086E+01, 0.5313E-02, 0.1265E+01, 0.8000E+01,
        0.2979E+00, 0.1354E+00, 0.2200E+00, 0.8296E-01, 0.9578E-01,
        0.1094E+00, 0.3794E+00, 0.8122E+01, 0.5189E+01, 0.3290E+01,
        0.1870E+01, 0.6110E+01,-0.3464E+02, 0.9000E+03, 0.1717E+01,
        0.4335E-01, 0.1915E+03, 0.2232E+00, 0.2119E+01, 0.2088E+01,
        -0.3029E+00, 0.2012E+00, 0.1104E-02, 0.2276E-01,-0.4562E+00,
        0.2397E+00, 0.1204E+01, 0.2321E-01, 0.5419E+03, 0.2247E+00,
        0.2168E+01, 0.2266E+03, 0.7649E-01, 0.1457E+01, 0.1318E+00,
        -0.7534E+02, 0.1776E+00, 0.1636E+01, 0.1350E+00,-0.5596E-02,
        0.5883E-02, 0.1934E+01, 0.3800E+00, 0.3319E+01, 0.1446E+00
    };

    double P[24] = {
        5.1377e-03,   9.8071e-01,   4.6379e-02,   1.6433e+00,
        6.9826e+00,  -2.2655e-01,   1.1095e-01,   2.7945e-02,
        4.0643e-01,   1.6076e+00,  -7.5460e+00,   4.4418e+00,
        -3.7464e-01,   1.0414e-01,  -2.6852e-01,   9.6653e-01,
        -1.9055e+00,   9.8965e-01,   2.0613e+02,  -4.5536e-02,
        2.4902e-01,  -1.3728e-01,   2.9201e+01,   4.9280e-03
    };

    double F1M;
    nu = (wsq - PM*PM + qsq) / 2. / PM;
    qv = sqrt(nu*nu + qsq);

    if (wsq <= 0.0) {
      W = 0.0;
      x = 0.0;
    }else{
      W  = sqrt(wsq);
      x  = qsq / (2.0 * PM * nu);
    }



//--------added by Yuxiang


    //cross section for proton or neutron
    W1=0;
    W2=0;

    if(IA<2 && wsq>1.155){
	    //cout<<"doing IA<2"<<endl;
	     christy507(wsq,qsq,F1p,Rc,sigt,sigl);
	     //if it is a neutron, subtract proton from deuteron, factor of two to 
	     //convert from per nucleon to per deuteron
	     if(Z<0.5){
		     resmodd(wsq,qsq, xvald0, F1d);
		     F1p = F1d*2.0 - F1p;
	     }

	     //cout<<"debug: "<<F1p<<"	"<<Rc<<endl;
	     W1=F1p/PM;
	     W2=W1*(1.0 + Rc)/(1.0+ nu*nu/qsq);
    }
//----------end added by Yuxiang
 


//------------------Yuxiang added deuteron----------------------
   double F1c=0;

   if(IA==2){
	   //get Fermi-smeared R from Erics proton fit
	   pind(wsq,qsq, F1c, Rc,sigt,sigl);
	   //get fit to F1 in deuteron, per nucleon
	   F1d=resd(qsq,wsq,xvald0);

	   W1=F1d/PM*2.0; //updated F1d

	   W2=W1*(1+Rc)/(1+nu*nu/qsq); // updated Rc
   }

//----------end Yuxiang added deuteron--------------------


//--------------
//For nuclei with A>2
    if(IA>2){

	    sigt = 0.;
	    sigl = 0.;
	    F1d = 0.;
	    F1p = 0.;

	    // Modifed to use Superscaling from Sick, Donnelly, Maieron,
	    // nucl-th/0109032
	    if(IA==2) {
	        kf=0.085;
	        Es=0.0022;
	    }else if(IA==3) {
	        kf=0.115;
	        Es=0.001;
	    }else if(IA>3 && IA<=7) {
	        kf=0.19;
	        Es=0.017;
	    }else if(IA>7 && IA<=16) {
	        kf=0.228;
	        Es=0.0165;
	    }else if(IA>16 && IA<=25) {
	        kf=0.230;
	        Es=0.025;
	    }else if(IA>25 && IA<=38) {
	        kf=0.236;
	        Es=0.018;
	    }else if(IA>38 && IA<=55) {
	        kf=0.241;
	        Es=0.028;
	    }else if(IA>55 && IA<=60) {
	        kf=0.241;
	        Es=0.018;
	    }else{ //at this point IA should be larger than 60   
		kf=0.245;
	        Es=0.018;
	    }
    
	    // adjust pf to give right width based on kf
	    pf = 0.5 * kf;
	    // assume this is 2 * pf * qv
	    DW2DPF = 2. * qv;
	    dw2des = 2. * (nu + PM) ;
	    // switched to using 99 bins!
	    for (ism = 0; ism < 99; ism++) {
	        Fyuse = fyp[ism]/100.;
	        Wsqp = wsq + XXp[ism] * pf * DW2DPF - Es * dw2des;
	        if (Wsqp > 1.159) {
	            christy507(Wsqp,qsq,F1pp,Rc,sigtp,siglp);
	            int bad=resmodd(Wsqp,qsq,xvald0,F1dp);
		    if(bad) return -1;
	            F1d = F1d + F1dp * Fyuse;
	            F1p = F1p + F1pp * Fyuse;
	            sigt = sigt + sigtp * Fyuse;
	            sigl = sigl + siglp * Fyuse;
	        }
	    }

	    Rc = 0.;
	    if(sigt > 0.) Rc = sigl / sigt;
	    W1 = (2. * Z * F1d + (IA - 2. * Z) * (2. * F1d - F1p)) / PM;
	    W1= W1*(1.0+P[13]*x+P[14]*pow(x,2)+P[15]*pow(x,3)+P[16]*pow(x,4)+P[17]*pow(x,5));
	    if(W > 0.0)  W1=W1*pow((1.0+(P[20]*W+P[21]*pow(W,2))/(1.0+P[22]*qsq)),2);
	    F1M = MEC2009(qsq,wsq,IA);
	    W1 = W1 + F1M;
	    if(wsq > 0.0 ) Rc = Rc * ( 1.0 + P[6] + P[23]*IA );
	    W2 = W1 * (1. + Rc) / (1. + nu*nu / qsq);
    } //end if(A>=2)

    x4 = qsq / 2. / PM / nu;
    emcfac = fitemc(x4, IA);
    F1 = PM * W1 * emcfac;
    F2 = nu * W2 * emcfac;
   /* if(std::isnan(F1) || std::isinf(F1) || std::isnan(F2) || std::isinf(F2)){
	    std::cout << "Error! nan/inf "<<__FILE__ << " line " << __LINE__ <<" "<<std::endl;
	    std::cout <<" "<< F1 <<" "<<F2<<" "<< W1 <<" "<<W2<<" "<<F1M<<" "<<PM<<" "<<nu<<" "<<emcfac<<std::endl;
	    return 1;
    }
    */
    return 0;
}








//-----------------------------------------------------------------//
//always calculates values
void christy507(double w2,double q2,double &F1,
                double &R, double &sigT, double &sigL)//checked
{
//   M.E. Christy and P.E. Bosted, ``Empirical Fit to Precision
//    Inclusive Electron-Proton Cross Sections in the Resonance Region'',
//    (arXiv:0712.3731). To be submitted to Phys. Rev. C.

    double xval1[50],xvalL[50];

    double mp = 0.938;
    double mp2 = mp*mp;
    double pi = 3.141593;
    double alpha = 1./137.036;

    double xval[100] = {
        0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
        0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
        0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
        0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
        0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
        -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
        0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
        0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
        -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
        0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
        0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
        0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
        0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
        0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
        0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
        0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
        0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
        0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
        0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
        0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00
    };

    for (int i = 0; i<50; i++) {
        xval1[i] = xval[i];
        xvalL[i] = xval[50+i] ;
        if(i < 12) xvalL[i] = xval1[i];
    }
    xvalL[42] = xval1[46];
    xvalL[43] = xval1[47];
    xvalL[49] = xval1[49];

    sigT = resmod507_v2(1,w2,q2,xval1);
    sigL = resmod507_v2(2,w2,q2,xvalL);

    F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3;
//     double xb = q2/(w2+q2-mp2);
    //double FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3;//seems to be not needed for anything
    R = sigL/sigT;

    return;
}








// -------------------------------------------------------------------------//
//on q2 or w2 out of range it returns -1
int resmodd(double w2, double q2,
	      double xval[50], double &sig) 
{
    //! returns F1 for average of free proton and neutron
    //! for given W2, Q2

    double W,mp,mp2,mass[7],width[7];
    double height[7],rescoef[6][4];
    double nr_coef[3][4],wdif,sig_nr;
    double sigr[7];
    double mpi,meta,intwidth[7],k,kcm,kr[7],kcmr[7],ppicm,ppi2cm;
    double petacm,ppicmr[7],ppi2cmr[7],petacmr[7],epicmr[7],epi2cmr[7];
    double eetacmr[7],epicm,epi2cm,eetacm,br[7][2],ang[7],sigrsv[7];
    double pgam[7],pwid[7][2],x0[7],dip;
    double sig_res,xpr,alpha,pi,F1;
    int i,j,num,iw;

    double xval0[12] = {
        0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
        0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
        0.10936E+00,0.37944E+00
    };

    double w2sv;

    mp = 0.938;
    mpi = 0.135;
    meta = 0.547;
    mp2 = mp*mp;
    pi = 3.141593;
    alpha = 1./137.036;

    sig = 0.;

    if(w2 < 1.07327*1.07327 || w2 > 25 || q2 < 0.0 || q2 > 11.0) {
	std::cout << "ERROR: " << __FILE__ << " line " << __LINE__ << std::endl;	
	std::cout << " W2/Q2 check failed: should be 1.07327^2<W2<25 and 0 < Q2 < 11" << std::endl;
	std::cout << "   but are w2 q2: " << w2 << " " << q2 << std::endl;
	return -1;
    }
    
    // ! branching ratios
    br[0][0] = 1.0;
    br[1][0] = 0.5;
    br[2][0] = 0.65;
    br[3][0] = 0.65;
    br[4][0] = 0.4;
    br[5][0] = 0.65;
    br[6][0] = 0.6;

    // ! angular momenta
    ang[0] = 1.;      // !!!  P33(1232)
    ang[1] = 0.;      // !!!  S11(1535)
    ang[2] = 2.;      // !!!  D13(1520)
    ang[3] = 3.;      // !!!  F15(1680)
    ang[4] = 0.;      // !!!  S15(1650)
    ang[5] = 1.;      // !!!  P11(1440) roper
    ang[6] = 3.;      // !!!  ? 4th resonance region

    // ! x0 parameter
    for (i = 0; i < 7; i++) x0[i] = 0.215;
    x0[0] = 0.1446;

    // ! out branching ratio
    for (i = 0; i < 7; i++) br[i][1] = 1.-br[i][0];

    // ! remember w2
    w2sv = w2;

    // ! uses xvals of 1-12, 47, and 48
    // ! move masses, wdiths into local variables
    // ! pyb changed to be fixed
    num = 0;
    for (i = 0; i < 6; i++) {
        num = num + 1;
        mass[i] = xval0[i];
    }
    for (i = 0; i < 6; i++) {
        num = num + 1;
        intwidth[i] = xval0[num-1];
    }
    
    //! changed to allow delta width, mass to vary
    //! taken out again since xval[1] used in MEC
    mass[6] = 1.9341;
    intwidth[6] = 0.380;

    iw = int(1000.*sqrt(w2));
    W = 0.001 * (iw+0.5);        
    w2 = W*W;
    wdif = W - (mp + mpi);
        
    // ! Calculate kinematics needed for threshold Relativistic B-W        
    k = (w2 - mp2) / 2. / mp;            
    kcm = (w2 - mp2) / 2. / W;    
    epicm = (w2 + mpi*mpi -mp2 ) / 2. / W;    
    ppicm = pow(std::max(0.0,(epicm*epicm - mpi*mpi)),0.5);    
    epi2cm = (w2 + pow((2.*mpi),2) -mp2 ) / 2. / W;    
    ppi2cm = pow(std::max(0.0,(epi2cm*epi2cm - pow((2.*mpi),2))),0.5);    
    eetacm = (w2 + meta*meta -mp2 ) / 2. / W;   
    petacm =  pow(std::max(0.0,(eetacm*eetacm - meta*meta)),0.5);            
    for (i = 0; i < 7; i++) {    
      kr[i] = (mass[i]*mass[i]-mp2)/2./mp;      
      kcmr[i] = (mass[i]*mass[i]-mp2)/2./mass[i];      
      epicmr[i] = (mass[i]*mass[i] + mpi*mpi -mp2 )/2./mass[i];      
      ppicmr[i] = pow(std::max(0.0,(epicmr[i]*epicmr[i] - mpi*mpi)),0.5);      
      epi2cmr[i] = (mass[i]*mass[i] + pow(2.*mpi,2) -mp2 )/2./mass[i];      
      ppi2cmr[i] = pow(std::max(0.0,(epi2cmr[i]*epi2cmr[i] - pow((2.*mpi),2))),0.5);      
      eetacmr[i] = (mass[i]*mass[i] + meta*meta -mp2 )/2./mass[i];      
      petacmr[i] =  pow(std::max(0.0,(eetacmr[i]*eetacmr[i] - meta*meta)),0.5);
            
      // ! Calculate partial widths      
      pwid[i][0] = intwidth[i]*pow((ppicm/ppicmr[i]),(2.*ang[i]+1.))*pow(((ppicmr[i]*ppicmr[i]+x0[i]*x0[i])/(ppicm*ppicm+x0[i]*x0[i])),ang[i]);            

      if(i != 1)       
	pwid[i][1] = intwidth[i]*pow((ppi2cm/ppi2cmr[i]),(2.*ang[i]+4.))*
		     pow(((ppi2cmr[i]*ppi2cmr[i]+x0[i]*x0[i])/(ppi2cm*ppi2cm+x0[i]*x0[i])),(ang[i]+2.))* W / mass[i];        	
      else  
	pwid[i][1] =  intwidth[1]*pow((petacm/petacmr[i]),(2.*ang[i]+1.))*
		      pow(((petacmr[i]*petacmr[i]+x0[i]*x0[i])/(petacm*petacm+x0[i]*x0[i])),ang[i]);        	
    
      pgam[i] = pow((kcm/kcmr[i]),2)*(kcmr[i]*kcmr[i]+x0[i]*x0[i])/(kcm*kcm+x0[i]*x0[i]);        	
      pgam[i] = intwidth[i]*pgam[i];        	
      width[i] = br[i][0]*pwid[i][0]+br[i][1]*pwid[i][1];        	
      sigr[i] = width[i] * pgam[i] / (pow(w2 - mass[i]*mass[i],2) + pow(mass[i]*width[i],2)) * kr[i] / k * kcmr[i] / kcm / intwidth[i];      
    }
    
    w2 = w2sv;

    // ! get parameters into local variables
    num = 12;
    //! resonance height coefficients. xvals of 13-36
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 4; j++) {
            num = num + 1;
            rescoef[i][j]=xval[num-1];
        }
    }
    // !  Non-Res coefficients xvals of 37-44
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 4; j++) {
            num = num + 1;
            nr_coef[i][j]=xval[num-1];
        }
    }

    // ! Begin resonance Q^2 dependence calculations   CCC
    // ! uses xvals 49
    for (i = 0; i < 6; i++) {
        height[i] = rescoef[i][0]*(1.+rescoef[i][1] * q2 / (1. + rescoef[i][2] * q2))/pow((1. + q2/0.91),rescoef[i][3]);
    }
    dip = 1./pow((1. + q2 / 0.91),2);
    height[6] = xval[48]*dip;
    sig_res = 0.;    
    for (i = 0; i < 7; i++) {
        sigrsv[i] =  height[i]*height[i] * sigr[i];
        sig_res = sig_res + sigrsv[i];
    }
    
    sig_res = sig_res * sqrt(w2);

    //! Begin non-resonant part uses xvals 45, 46, 50
    //! Depends on both W2 and Q2 so can't easily precalculate
    sig_nr = 0.;
    xpr = 1.+(w2-pow(mp+mpi,2))/(q2+0.05);
    xpr = 1./xpr;
    W = pow(w2,0.5);
    wdif = W - (mp + mpi);
    for (i = 0; i < 2; i++) {
        sig_nr = sig_nr +(nr_coef[i][0]*pow(wdif,((2*(i+1)+1)/2.)))/pow(q2+nr_coef[i][1],(nr_coef[i][2]+nr_coef[i][3]*q2+xval[44+i]*q2*q2));
    }
    sig_nr = sig_nr * xpr;
    sig = sig_res + sig_nr;

    F1 = sig * (w2-mp2)/8./pi/pi/alpha/0.3894e3;
    sig = F1;

    return 0;
}










//------------------------------------------------------------------//
// Used for Aluminum inelastic...
// always returns a calculated value (called resmod507 in the fortran code)
double resmod507_v2(double sf,double w2,
                      double q2,double xval[50])
{

    double W,mp,mp2,xb;
    double mass[7] = {0,0,0,0,0,0,0};
    double width[7] = {0,0,0,0,0,0,0};
    double height[7] = {0,0,0,0,0,0,0};
    double intwidth[7] = {0,0,0,0,0,0,0};
    double rescoef[6][4];
    double nr_coef[3][4],sigr[7],wdif[2],sig_nr,h_nr[3];
    double mpi,meta,k,kcm,kr[7],kcmr[7],ppicm,ppi2cm;
    double petacm,ppicmr[7],ppi2cmr[7],petacmr[7],epicmr[7],epi2cmr[7];
    double eetacmr[7],epicm,epi2cm,eetacm,br[7][3],ang[7];
    double pgam[7],pwid[7][3],x0[7],q20;
    double sig_res,t,xpr[2],m0,sig;
    int i,j,num;

    mp = 0.938;
    mpi = 0.135;
    meta = 0.547;
    mp2 = mp*mp;
    W = sqrt(w2);
    wdif[0] = W - (mp + mpi);
    wdif[1] = W - (mp + 2.*mpi);

    m0 = 0.125;
    if(sf == 2) m0 = xval[48];
    if(sf == 1) {
        q20 = 0.05;
    } else {
        q20 = 0.125;
    }

    // single pion branching ratios  CCCC

    br[0][0] = 1.0;       //!!!  P33(1232)
    br[1][0] = 0.45;      //!!!  S11(1535)
    br[2][0] = 0.65;      //!!!  D13(1520)
    br[3][0] = 0.65;      //!!!  F15(1680)
    br[4][0] = 0.4;       //!!!  S11(1650)
    br[5][0] = 0.65;      //!!!  P11(1440) roper
    br[6][0] = 0.50 ;     //!!!  F37(1950)

    br[0][2] = 0.0;       //!!!  P33(1232)
    br[1][2] = 0.45;      //!!!  S11(1535)
    br[2][2] = 0.0;       //!!!  D13(1520)
    br[3][2] = 0.0;       //!!!  F15(1680)
    br[4][2] = 0.1;       //!!!  S11(1650)
    br[5][2] = 0.0;       //!!!  P11(1440) roper
    br[6][2] = 0.0;       //!!!  F37(1950)

    // 2-pion branching ratios  CCCC

    for (i = 0; i < 7; i++) {
        br[i][1] = 1.-br[i][0]-br[i][2];
    }

    // Meson angular momentum   CCCC
    ang[0] = 1.;       //!!!  P33(1232)
    ang[1] = 0.;       //!!!  S11(1535)
    ang[2] = 2.;       //!!!  D13(1520)
    ang[3] = 3.;       //!!!  F15(1680)
    ang[4] = 0.;       //!!!  S15(1650)
    ang[5] = 1.;       //!!!  P11(1440) roper
    ang[6] = 3.;       //!!!  F37(1950)

    for (i = 0; i < 7; i++) { //!!!  resonance damping parameter  !!!
        x0[i] = 0.215;
    }
    x0[0] = 0.15;
    x0[0] = xval[49];

    xb = q2/(q2+w2-mp2);
    xpr[0] = 1.+(w2-pow((mp+mpi),2))/(q2+q20);
    xpr[0] = 1./xpr[0];
    xpr[1] = 1.+(w2-pow((mp+mpi+mpi),2))/(q2+q20);
    xpr[1] = 1./xpr[1];

    t = log(log((q2+m0)/0.330/0.330)/log(m0/0.330/0.330));

    // Calculate kinematics needed for threshold Relativistic B-W

    k = (w2 - mp2)/2./mp;
    kcm = (w2-mp2)/2./W;

    epicm = (w2 + mpi*mpi -mp2 )/2./W;
    ppicm = pow(std::max(0.0,(epicm*epicm - mpi*mpi)),0.5);
    epi2cm = (w2 + pow((2.*mpi),2) -mp2 )/2./W;
    ppi2cm = pow(std::max(0.0,(epi2cm*epi2cm - pow((2.*mpi),2))),0.5);
    eetacm = (w2 + meta*meta -mp2 )/2./W;
    petacm =  pow(std::max(0.0,(eetacm*eetacm - meta*meta)),0.5);

    num = 0;

    for (i = 0; i < 6; i++) { // !!!  Read in resonance masses     !!!
        num = num + 1;
        mass[i] = xval[i];
    }
    for (i = 0; i < 6; i++) { // !!!  Read in resonance widths     !!!
        num = num + 1;
        intwidth[i] = xval[num-1];
        width[i] = intwidth[i];
    }
    if (sf == 2) { //      !!!  Put in 4th resonance region  !!!
        mass[6] = xval[42];
        intwidth[6] = xval[43];
        width[6] = intwidth[6];
    } else {
        mass[6] = xval[46];
        intwidth[6] = xval[47];
        width[6] = intwidth[6];
    }

    for (i = 0; i < 7; i++) {
        kr[i] = (mass[i]*mass[i]-mp2)/2./mp;
        kcmr[i] = (mass[i]*mass[i]-mp2)/2./mass[i];
        epicmr[i] = (mass[i]*mass[i] + mpi*mpi -mp2 )/2./mass[i];
        ppicmr[i] = pow(std::max(0.0,(epicmr[i]*epicmr[i] - mpi*mpi)),0.5);
        epi2cmr[i] = (mass[i]*mass[i] + pow((2.*mpi),2) -mp2 )/2./mass[i];
        ppi2cmr[i] = pow(std::max(0.0,(epi2cmr[i]*epi2cmr[i] - pow((2.*mpi),2))),0.5);
        eetacmr[i] = (mass[i]*mass[i] + meta*meta -mp2 )/2./mass[i];
        petacmr[i] =  pow(std::max(0.0,(eetacmr[i]*eetacmr[i] - meta*meta)),0.5);

        // CCC   Calculate partial widths   CCC

        pwid[i][0] = intwidth[i]*pow((ppicm/ppicmr[i]),(2.*ang[i]+1.))*pow(((ppicmr[i]*ppicmr[i]+x0[i]*x0[i])/(ppicm*ppicm+x0[i]*x0[i])),ang[i]); //      !!!  1-pion decay mode

        pwid[i][1] = intwidth[i]*pow((ppi2cm/ppi2cmr[i]),(2.*ang[i]+4.))*pow(((ppi2cmr[i]*ppi2cmr[i]+x0[i]*x0[i])/(ppi2cm*ppi2cm+x0[i]*x0[i])),(ang[i]+2)); //  !!!  2-pion decay mode

        pwid[i][1] = W/mass[i]*pwid[i][1];
        pwid[i][2] = 0.; //          !!!  eta decay mode

        if(i == 1 || i == 4) {
            pwid[i][2] =  intwidth[i]*pow((petacm/petacmr[i]),(2.*ang[i]+1.))*pow(((petacmr[i]*petacmr[i]+x0[i]*x0[i])/(petacm*petacm+x0[i]*x0[i])),ang[i]); // !!!  eta decay only for S11's
        }

        pgam[i] = pow((kcm/kcmr[i]),2)*(pow(kcmr[i],2)+x0[i]*x0[i])/(kcm*kcm+x0[i]*x0[i]);
        pgam[i] = intwidth[i]*pgam[i];
        width[i] = br[i][0]*pwid[i][0]+br[i][1]*pwid[i][1]+br[i][2]*pwid[i][2];
    }

    //CCC    End resonance kinematics and Widths calculations   CCC

    // CCC    Begin resonance Q^2 dependence calculations   CCC

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 4; j++) {
            num = num + 1;
            rescoef[i][j]=xval[num-1];
        }

        if(sf == 1) {
            height[i] = rescoef[i][0]*(1.+rescoef[i][1]*q2/(1.+rescoef[i][2]*q2))/pow((1.+q2/0.91),rescoef[i][3]);
        } else {
            height[i] = rescoef[i][0]*q2/(1.+rescoef[i][1]*q2)*exp(-1.*rescoef[i][2]*q2);
        }

        height[i] = height[i]*height[i];
    }

    if(sf == 2) { // !!!  4th resonance region  !!!
        height[6] = xval[44]*q2/(1.+xval[45]*q2)*exp(-1.*xval[46]*q2);
    } else {
        height[6] = xval[48]/(1.+q2/0.91);
    }

    height[6] = height[6]*height[6];

    // CCC    End resonance Q^2 dependence calculations   CCC

    for (i = 0; i < 3; i++) { // !!!  Non-Res coefficients  !!!
        for (j = 0; j < 4; j++) {
            num = num + 1;
            nr_coef[i][j]=xval[num-1];
        }
    }

    // CCC   Calculate Breit-Wigners for all resonances   CCC

    sig_res = 0.0;
    for (i = 0; i < 7; i++) {
        sigr[i] = width[i]*pgam[i]/(pow(w2 - mass[i]*mass[i],2.) + pow(mass[i]*width[i],2.));
        sigr[i] = height[i]*kr[i]/k*kcmr[i]/kcm*sigr[i]/intwidth[i];
        sig_res = sig_res + sigr[i];
    }
    sig_res = sig_res*W;
    // CCC    Finish resonances / start non-res background calculation   CCC
    sig_nr = 0.;
    if(sf == 1) {
        for (i = 0; i < 2; i++) {
            h_nr[i] = nr_coef[i][0]/pow(q2+nr_coef[i][1],nr_coef[i][2] + nr_coef[i][3]*q2+xval[44+i]*q2*q2);
            sig_nr = sig_nr +h_nr[i]*(pow(wdif[0],(2*(i+1)+1)/2.));
        }

        sig_nr = sig_nr*xpr[0];

    } else if (sf == 2) {
        for (i = 0; i < 1; i++) {
            sig_nr = sig_nr + nr_coef[i][0]*pow((1.-xpr[i]),(nr_coef[i][2]+nr_coef[i][1]*t))/(1.-xb)/(q2+q20);
            sig_nr = sig_nr*pow(q2/(q2+q20),nr_coef[i][3])*pow(xpr[i],(xval[40]+xval[41]*t));
        }
    }
    sig = sig_res + sig_nr;

    return sig;
}












double MEC2009(double q2,double w2, int A)//checked
{
    //! fit to low q2 dip region: purefly empirical
    //! assume contribution is pure transverse
    //returns 0 on A<2.5 and W2<=0, otherwise returns calculated value
    //propagated to F1F2IN09 into W1 as W1+=MEC2009
    double f1 = 0.0;
    double am = 0.938;
    double w,nu;

    double P[24] = {
        5.1377e-03,   9.8071e-01,   4.6379e-02,   1.6433e+00,
        6.9826e+00,  -2.2655e-01,   1.1095e-01,   2.7945e-02,
        4.0643e-01,   1.6076e+00,  -7.5460e+00,   4.4418e+00,
        -3.7464e-01,   1.0414e-01,  -2.6852e-01,   9.6653e-01,
        -1.9055e+00,   9.8965e-01,   2.0613e+02,  -4.5536e-02,
        2.4902e-01,  -1.3728e-01,   2.9201e+01,   4.9280e-03
    };

    double p18, x, f1corr;

    if (w2 <= 0.0) return 0;
    w  = sqrt(w2);
    nu = (w2 - am*am + q2) / 2. / am;
    x  = q2 / (2.0 * am * nu );

    if (A < 2.5) return 0;

    p18 = P[18];
    // ! special case for 3He
    if (A > 2.5 && A < 3.5) p18 = 70;
    // ! special case for 4He
    if (A > 3.5 && A < 4.5) p18 = 170;
    // ! new values for C, Al, Cu
    if(A > 4.5) p18 = 215;
    if(A > 20.) p18 = 235;
    if(A > 50.) p18 = 230;


    f1corr = P[0]*exp(-pow(w-P[1],2)/P[2])/(pow(1.0 + std::max(0.3,q2)/P[3],P[4]))*pow(nu,P[5])*(1.0 + p18 * pow(A,(1.0 + P[19] * x)));

    f1 = f1corr;

    if (f1 <= 1.0E-9 ) f1 = 0.0;

   /* if(std::isnan(f1)){
      std::cout << "Error! nans "<<__FILE__ << " line " << __LINE__ <<" "<<std::endl;
      std::cout <<" "<< f1 <<" "<<x<<" "<< nu <<" "<<w2<<" "<<p18<<" "<<q2<<" "<<std::endl;
      std::cout <<"  "<< exp(-pow(w-P[1],2)/P[2]) << " " << (pow(1.0 + std::max(0.3,q2)/P[3],P[4])) <<" "<< pow(nu,P[5]) 
	     << " " << (1.0 + p18 * pow(A,(1.0 + P[19] * x)))<< std::endl;
    }
    */
    return f1;
}






double fitemc(double X, int A)
{
    /*!---------------------------------------------------------------------
      ! Fit to EMC effect.  Steve Rock 8/3/94
      ! Funciton returns value of sigma(A)/sigma(d)
      ! with no isoscalerity correction
      ! A= atomic number
      ! x = Bjorken x.
      !
      ! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number
      ! First data at each x was fit to form C*A**alpha.  The point A=2(d)
      !  was includded with a value of 1 and an error of about 2%.
      ! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.
      ! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91)
      !  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
      ! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
      ! C(x) was fit to a 3 term polynomial in natural logs as
      !  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.

      ! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
      !                    also gave value at x=.0085 if x>.88
      ! 11/05 PYB PEB modified to use value at x=0.7 if x>0.7, because
      !    beyond that assume Fermi motion is giving the rise, and we
      !    already are taking that into account with the y-smearing of
      !    the inelastic
      !-----------------------------------------------------------------------
    */
  // according to fortran 2009 code: return 1 (no emc effects) on A<2.5 otherwise returns calculated value
  //propagated into F1F2IN09 for the calculation of F1 as a multiplicative factor

    double ALPHA,C,LN_C,X_U;
    double fitval;

    // !Chisq=         19.   for 30 points
    // !Term    Coeficient
    double ALPHA_COEF[9] = {
        -6.98871401E-02,
        2.18888887E+00,
        -2.46673765E+01,
        1.45290967E+02,
        -4.97236711E+02,
        1.01312929E+03,
        -1.20839250E+03,
        7.75766802E+02,
        -2.05872410E+02
    };

    // !Chisq=         22.    for 30 points
    // !Term    Coeficient
    double C_COEF[3] = {
        1.69029097E-02,
        1.80889367E-02,
        5.04268396E-03
    };

    fitval = 1.;
    if (A < 2.5) return 1;

    if ( (X > 0.70) || (X < 0.0085) ) {
        // !Out of range of fit
        if (X < 0.0085) X_U =.0085;
        if (X > 0.70) X_U = 0.70;
    } else {
        X_U = X;
    }

    LN_C = C_COEF[0];
    for (int i = 1; i <= 2; i++) {
        LN_C = LN_C + C_COEF[i] * pow(log(X_U),i);
    }

    C = exp(LN_C);

    ALPHA = ALPHA_COEF[0];
    for (int i = 1; i <= 8; i++) {
        ALPHA = ALPHA + ALPHA_COEF[i] * pow(X_U,i);
    }

    fitval = C*pow(A,ALPHA);
    return fitval;
}







//#####################################################################################
double pind(double w2, double q2, double &F1, double &R, double &sigt, double &sigl){
	// calculate proton with Fermi smearing fo a deuteron
	
	double am=0.938;
	double nu, qv, F1p, Rp, sigtp, siglp;
	double amd=1.8756;
	double w2p, pz;
       
       double fyd[20]={0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041,
     		       0.5029, 0.5034, 0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992,
     		       0.4994, 0.4977, 0.5023, 0.4964, 0.4966, 0.4767};

       double avpz[20]={-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,
     			 -0.0195,-0.0135,-0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199,
      			 0.0268, 0.0349, 0.0453, 0.0598, 0.0844, 0.1853};
       
       double avp2[20]= {0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068,
     			 0.0060, 0.0054, 0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060,
     			 0.0069, 0.0081, 0.0102, 0.0140, 0.0225, 0.0964};

	// Look up tables for deuteron in fine bins for sub threshold
       double fydf[200]={
      0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
      0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
      0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
      0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
      0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
      0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
      0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
      0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
      0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
      0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
      0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
      0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
      5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
      4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
      0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
      0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
      0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
      0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
      0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
      0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
      0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
      0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
      0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
      0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
      0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001};

       double avp2f[200]={
          1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
      0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
      0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
      0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
      0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
      0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
      0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
      0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
      0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
      0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
      0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
      0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
      0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
      0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
      0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
      0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
      0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
      0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
      0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
      0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
      0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
      0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
      0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
      0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
      0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0};

      nu = (w2 - am*am + q2) / 2. / am;
      qv = sqrt(nu*nu + q2);
      F1=0.;
      R=0.;
      sigt=0.;
      sigl=0.;
      
      //Do fast 20 bins if abvoe threshold
      if(w2>1.3){
	      for(int ism=0; ism<20; ism++){
		      w2p=(amd + nu - sqrt(am*am + avp2[ism])) * (amd + nu - sqrt(am*am + avp2[ism])) - qv*qv + 2.0*qv*avpz[ism] - avp2[ism];
		      if(w2p>1.155){
			      christy507(w2p, q2, F1p, Rp, sigtp, siglp);
			      sigt=sigt + sigtp*fyd[ism]/10.0;
			      sigl=sigl + siglp*fyd[ism]/10.0;
			      F1 = F1 + F1p*fyd[ism]/10.0;
		      }
	      }
      }else{
	      for(int ism=0; ism<200; ism++){
		      pz=-1.0 + 0.01*(ism-0.5);
		      w2p = (amd +nu - sqrt(am*am + avp2f[ism]) ) * (amd +nu - sqrt(am*am + avp2f[ism]) ) - qv*qv + 2*qv*pz-avp2f[ism];
		      if(w2p>1.155){
			      christy507(w2p, q2, F1p, Rp,sigtp, siglp);
			      sigt = sigt+sigtp*fydf[ism]/100.0;
			      sigl = sigl + siglp*fydf[ism]/100.0;
			      F1 = F1 + F1p*fydf[ism]/100.0;
		      }
	      }
      }

      if(sigt!=0)R=sigl/sigt;

      return 0;
}








//#####################################################################
double resd(double q2, double w2, double xval[50]){
	//return F1

	//calculate deuteron F1 by Fermi smearing of proton plus neutron
	//add MEC term (not smeared)
	//arxiv 0711.0159, submitted to PRC
	//epirical fit to inelastic electron-deuteron and electron-neutron resonance 
	//region transverse cross sections.
	
	double F1, nu, qv, dw2dpf, w2p, sigp, f1sv;
	double am=0.938;
	double amd=1.8756, pz, f1m, f2m;

        double fyd[20]= {0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041,
       0.5029, 0.5034, 0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992,
       0.4994, 0.4977, 0.5023, 0.4964, 0.4966, 0.4767};

        double avpz[20]={-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,
      -0.0195,-0.0135,-0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199,
       0.0268, 0.0349, 0.0453, 0.0598, 0.0844, 0.1853};

       double avp2[20]={ 0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068,
       0.0060, 0.0054, 0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060,
       0.0069, 0.0081, 0.0102, 0.0140, 0.0225, 0.0964};

       // Look up tables for deuteron in fine bins for sub threshold
       double fydf[200]={
      0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
      0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
      0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
      0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
      0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
      0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
      0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
      0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
      0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
      0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
      0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
      0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
      5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
      4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
      0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
      0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
      0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
      0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
      0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
      0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
      0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
      0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
      0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
      0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
      0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001};

      double avp2f[200]={
          1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
      0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
      0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
      0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
      0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
      0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
      0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
      0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
      0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
      0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
      0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
      0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
      0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
      0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
      0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
      0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
      0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
      0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
      0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
      0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
      0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
      0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
      0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
      0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
      0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0};


      nu = (w2 - am*am + q2) / 2. / am;
      qv = sqrt(nu*nu + q2);
      F1 = 0.;
      // Do fast 20 bins if abvoe threshold
      if(w2>1.3){
	      for(int ism=0; ism<20;ism++){
		      w2p=( amd + nu -sqrt(am*am + avp2[ism]) ) * ( amd + nu -sqrt(am*am + avp2[ism]) ) - qv*qv + 2.0*qv*avpz[ism] - avp2[ism];
		      if(w2p>1.155){
			      resmodd(w2p, q2, xval, sigp);
			      F1=F1 + sigp*fyd[ism]/10.0;
		      }
	      }
      }else{
	      F1=0;
	      for(int ism=0; ism<200; ism++){
		      pz=-1 + 0.01*(ism-0.5);
		      w2p=(amd + nu - sqrt(am*am + avp2f[ism]) ) * (amd + nu - sqrt(am*am + avp2f[ism]) ) - qv*qv + 2.0*qv*pz - avp2f[ism];
		      if(w2p>1.155){
			      resmodd(w2p,q2,xval,sigp);
			      F1=F1+sigp*fydf[ism]/100.0;
		      }
	      }
      }

      //add MEC term
      //took out again, then put back in again
      mec(1,2,q2,w2,f1m,f2m,xval);
      F1=F1 + f1m;

      return F1;
}





double mec(double z, double a, double q2, double w2, double &f1, double &f2, double xval[50]){
	double am=0.938;
	double w;
	double nu;

	w=sqrt(w2);
	nu = (w2 - am*am +q2)/2.0/am;
	
	double max_q2=0;
	if(q2>=0.3) max_q2=q2;
	if(q2<0.3) max_q2=0.3;
	
	f1= xval[0]*exp(-(w-xval[1])*(w-xval[1])/xval[2] )/pow( 1+ max_q2/xval[3], xval[4]) * pow(nu, xval[5]);

	f2=0;
	if(q2>0) f2=nu*(f1/am)/(1+nu*nu/q2);

	return 0;
}


//oooooooooooooooooooooo000ooooooooooooooooooooooooo
// Qweak Quasi Elastic part of the code
//oooooooooooooooooooooo000ooooooooooooooooooooooooo

void F1F2QE09(int Z, int IA, double QSQ,
	      double wsq, double &F1, double &F2)
{
//=======================================================================
//      SUBROUTINE F1F2QE09(Z, A, QSQ, wsq, F1, F2)
//
// Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS
// for A>2 uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
// for A=2 uses pre-integrated Paris wave function (see ~bosted/smear.f)
// coded by P. Bosted August to October, 2006
//
// input: Z, A  (real*8) Z and A of nucleus (shoud be 2.0D0 for deueron)
//        Qsq (real*8) is 4-vector momentum transfer squared (positive in
//                     chosen metric)
//        Wsq (real*8) is invarinat mass squared of final state calculated
//                     assuming electron scattered from a free proton
//                 
// outputs: F1, F2 (real*8) are structure functions per nucleus
//
// Note: Deuteron agrees well with Laget (see ~bosted/eg1b/laget.f) for
// a) Q2<1 gev**2 and dsig > 1% of peak: doesnt describe tail at high W
// b) Q2>1 gev**2 on wings of q.e. peak. But, this model is up
//    to 50% too big at top of q.e. peak. BUT, F2 DOES agree very
//    nicely with Osipenko et al data from CLAS, up to 5 GeV**2
  
//  cout << "Z, A, Q2, W2: " << Z << ", " << IA << ", " << QSQ << ", " << wsq << endl;

  double avgN, Pauli_sup1, Pauli_sup2, GEP, GEN, GMP, GMN, Q, Q3, Q4;
  double amp = 0.938;
  double amd=1.8756;
  double RMUP = 2.792782;
  double RMUN = -1.913148;
  double Nu;
  double QV, TAU, FY;
  double kappa, lam, lamp, taup, squigglef, psi, psip, nuL, nuT;
  double kf, Es, GM2bar, GE2bar, W2bar, Delta, GL, GT;
	
  int izz, izzmin, izznom, izdif;
  double pz, dpz, pznom, pzmin, dwin, w2p;
  
  // Look up tables for deuteron case
       double  fyd[200]={
      0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
      0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
      0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
      0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
      0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
      0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
      0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
      0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
      0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
      0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
      0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
      0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
      5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
      4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
      0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
      0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
      0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
      0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
      0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
      0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
      0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
      0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
      0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
      0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
      0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001};

       double avp2[200]={
          1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
      0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
      0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
      0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
      0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
      0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
      0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
      0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
      0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
      0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
      0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
      0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
      0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
      0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
      0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
      0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
      0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
      0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
      0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
      0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
      0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
      0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
      0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
      0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
      0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0};

  
  // Peter Bosted's correction params
  /*
  double pb[20] = {
    0.1023E+02, 0.1052E+01, 0.2485E-01, 0.1455E+01,
    0.5650E+01,-0.2889E+00, 0.4943E-01,-0.8183E-01,
    -0.7495E+00, 0.8426E+00,-0.2829E+01, 0.1607E+01,
    0.1733E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
    0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00};
  */
  double y,R;

  double P[24] = {
     5.1377e-03,   9.8071e-01,   4.6379e-02,   1.6433e+00,
     6.9826e+00,  -2.2655e-01,   1.1095e-01,   2.7945e-02,
     4.0643e-01,   1.6076e+00,  -7.5460e+00,   4.4418e+00,
     -3.7464e-01,   1.0414e-01,  -2.6852e-01,   9.6653e-01,
     -1.9055e+00,   9.8965e-01,   2.0613e+02,  -4.5536e-02,
     2.4902e-01,  -1.3728e-01,   2.9201e+01,   4.9280e-03};

  // return if proton: future change this to allow for
  // equivalent W resolution
  F1 = 0.;
  F2 = 0.;
  avgN = IA - Z;
  if (IA==1) return;

  // some kinematic factors. Return if Nu or QSQ is negative
  Nu = (wsq - amp*amp + QSQ) / 2. / amp;

  //cout << "In call... IA, Nu, QSQ = " << IA << ", " << Nu << ", " << QSQ << endl;
  if(Nu <= 0.0 || QSQ < 0.) return;
  TAU   = QSQ / 4.0 / amp / amp;                                 
  QV = sqrt(Nu*Nu + QSQ);

  // Bosted fit for nucleon form factors Phys. Rev. C 51, p. 409 (1995)
  Q = sqrt(QSQ);
  Q3 = QSQ * Q;
  Q4 = QSQ*QSQ;
  GEP = 1./  (1. + 0.14 * Q + 3.01 * QSQ + 0.02 * Q3 + 1.20 * Q4 + 0.32 * pow(Q,5));
  GMP = RMUP * GEP;
  GMN = RMUN / (1.- 1.74 * Q + 9.29 * QSQ - 7.63 * Q3 + 4.63 * Q4);
  GEN = 1.25 * RMUN * TAU / (1. + 18.3 * TAU) / pow((1. + QSQ / 0.71),2);

  //cout << "Form Factors: " << GEP << ", " << GMP << ", " << GEN << ", " << GMN << endl;

  // Get kf and Es from superscaling from Sick, Donnelly, Maieron,
    // nucl-th/0109032
    if(IA==2) kf=0.085;
    if(IA==2) Es=0.0022;
    // changed 4/09
    if(IA==3) kf=0.115;
    if(IA==3) Es=0.001 ;
    // changed 4/09
    if(IA>3) kf=0.19;
    if(IA>3) Es=0.017; 
    if(IA>7) kf=0.228;
    if(IA>7) Es=0.020;
    // changed 5/09
    if(IA>7) Es=0.0165;
    if(IA>16) kf=0.230;
    if(IA>16) Es=0.025;
    if(IA>25) kf=0.236;
    if(IA>25) Es=0.018;
    if(IA>38) kf=0.241;
    if(IA>38) Es=0.028;
    if(IA>55) kf=0.241;
    if(IA>55) Es=0.023;
    if(IA>60) kf=0.245;
    if(IA>60) Es=0.028;
    // changed 5/09 
    if(IA>55) Es=0.018;

  // Pauli suppression model from Tsai RMP 46,816(74) eq.B54
  if ((QV > 2.* kf) || (IA == 1)) {
    Pauli_sup2 =1.0;
  } else {
    Pauli_sup2 = 0.75 * (QV / kf) * (1.0 - (pow((QV / kf),2))/12.);
  }
  Pauli_sup1 = Pauli_sup2;

  //cout << "kf, Es, Paulisup1,2: " << kf << ", " << Es << ", " << Pauli_sup1 << ", " << Pauli_sup2 << endl;
  
  // structure functions with off shell factors
  kappa = QV / 2. / amp;
  lam = Nu / 2. / amp;
  lamp = lam - Es / 2. / amp;
  taup = kappa*kappa - lamp*lamp;
  squigglef = sqrt(1. + pow((kf/amp),2)) -1.;

  // Very close to treshold, could have a problem
  if(1.+lamp <= 0.) return;
  if(taup * (1. + taup) <= 0.) return;

  psi =  (lam  - TAU ) / sqrt(squigglef) / sqrt((1.+lam )* TAU + kappa * sqrt(TAU * (1. + TAU)));
  psip = (lamp - taup) / sqrt(squigglef) / sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)));
  nuL = pow((TAU / kappa / kappa),2);

  // changed definition of nuT from
  // nuT = TAU / 2. / kappa**2 + tan(thr/2.)**2
  // to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
  nuT = TAU / 2. / kappa / kappa;

  GM2bar = Pauli_sup1 * (Z * GMP*GMP + avgN * GMN*GMN);
  GE2bar = Pauli_sup2 * (Z * GEP*GEP + avgN * GEN*GEN);
  //double W1bar = TAU * GM2bar;
  W2bar = (GE2bar + TAU * GM2bar) / (1. + TAU);

  Delta = squigglef * (1. - psi*psi) * (sqrt(TAU * (1.+TAU)) / kappa + squigglef/3. *
     (1. - psi*psi) * TAU / kappa / kappa);
  GL = kappa*kappa / TAU * (GE2bar + Delta * W2bar) / 2. / kappa / (1. + squigglef * 
     (1. + psi*psi) / 2.);
  GT = (2. * TAU * GM2bar + Delta * W2bar) / 2. / kappa / (1. + squigglef * 
     (1. + psi*psi) / 2.);

  //cout << "nuL, nuT, GL, GT: " << nuL << ", " << nuT << ", " << GL << ", " << GT << endl;

  // added to prevent negative xsections:
  if (GT < 0) {
    GT = 0;
    //cout << "Reset GT to zero" << endl;
  }

  // from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
  FY = 1.5576 / (1. + 1.7720*1.7720 * pow((psip + 0.3014),2)) / (1. + exp(-2.4291 * psip)) / kf;




// USE PWIA and Paris W. F. for deuteron to get better FY 

  if(IA == 2){
	  pz=(QSQ - 2*amp*Nu)/2.0/QV;
	  izz = int((pz+1.0)/0.01)+1;
	  int max_izz;
	  if(izz>1){
		  max_izz=izz;
	  }else if(izz<=1){
		  max_izz=1;
	  }

	  if(max_izz>200){
		  izz=200;
	  }else if(max_izz<=200){
		  izz=izz;
	  }

	  izznom=izz;

	  //ignoring energy term, estimate change in pz to compensate
	  //for avp2 term
	  dpz=avp2[izznom-1]/2.0/QV;
	  izdif=dpz*150.0;
	  dwin=1.0e6;
	  izzmin=0;
	  
	  int loop_upper_max;
	  if((izznom+izdif)>1){
		  loop_upper_max=izznom+izdif;
	  }else if((izznom+izdif)<=1){
		  loop_upper_max=1;
	  }

	  int loop_upper;
	  if(loop_upper_max>200){
		  loop_upper=200;
	  }else if(loop_upper_max<=200){
		  loop_upper=loop_upper_max;
	  }

	  for(int izp=izznom;izp<=loop_upper;izp++){
		  pz=-1.0 + 0.01*(izp-0.5);
		  w2p=(amd+Nu-amp)*(amd+Nu-amp) - QV*QV + 2*QV*pz - avp2[izp-1];
		  if(abs(w2p-amp*amp) > dwin){
			  break;
		  }else if(abs(w2p-amp*amp)<=dwin){
			  dwin=abs(w2p-amp*amp);
			  izzmin=izp;
		  }
	  }

	  int tmp;
	  if(izzmin>2){
		tmp=izzmin;
	  }else if(izzmin<=2){
		tmp=2;
	  }

	  if(tmp>199){
		izz=199;
	  }else if(tmp<=199){
		izz=tmp;
	  }
			
	  pznom=-1 + 0.01*(izz-0.5);
	  dwin=1.e6;
	  for(int izp=1;izp<=19;izp++){
		  pz=pznom-0.01+0.001*izp;
		  w2p=(amd+Nu-amp)*(amd+Nu-amp)-QV*QV + 2*QV*pz- avp2[izz-1];
		  if(abs(w2p - amp*amp)<dwin){
			  dwin=abs(w2p-amp*amp);
			  pzmin=pz;
		  }
	  }

	  if(pzmin<pznom){
		  FY=fyd[izz-1] - (fyd[izz-2] - fyd[izz-1])*(pzmin-pznom)/0.01;
	  }else{
		  FY=fyd[izz-1]+ (fyd[izz+1-1]- fyd[izz-1])*(pzmin-pznom)/0.01;
	  }
  }

	  



  // final results
  F2 = Nu * FY * (nuL * GL + nuT * GT);
  F1 = amp * FY * GT / 2.;

  //cout << "nu, Fy, nuL, GL, nuT, GT, amp: " << endl;
  //cout << Nu << ", " << FY << ", " << nuL << ", " << GL << ", " << nuT << ", " << GT << ", " << amp << endl;

  if (F1 < 0.0) F1 = 0.;
  if (Nu > 0. && F1 > 0.) R = (F2 / Nu) / (F1 / amp) * (1. + Nu*Nu / QSQ) - 1.0;
  else R = 0.4/QSQ;


  // apply correction factors
  if ( IA > 2 ) {
    y = (wsq -amp*amp) / QV;
    //         F1 = F1 * (1. + pb(8) + pb(9) * y +
    //     >        pb(10)*y**2 + pb(11)*y**3 + pb(12)*y**4 )
    //         R = R * (1. + pb(13))
    //         F2 = Nu * F1/amp * (1. + R) / (1. + Nu**2/QSQ)

    // correction to correction Vahe
    if (wsq > 0.0) {
      F1=F1*(1.0+P[7]+P[8]*y+P[9]*y*y +P[10]*pow(y,3) +P[11]*pow(y,4));
      R = R * ( 1.0 + P[12] );
      F2 = Nu * F1/amp * (1. + R) / (1. + Nu*Nu/QSQ);
      if (F1 < 0.0) F1=0.0;
    }
  }
}










