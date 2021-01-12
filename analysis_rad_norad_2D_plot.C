#include <iostream> 
#include <fstream>
#include <cmath> 
#include "math.h" 
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TText.h"
#include "TSystem.h"
#include "TArc.h"
#include "TString.h"
#include <vector>
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"
using namespace std;
//#####################################################################################################################################################

int analysis_rad_norad_2D_plot(){
gStyle->SetPaintTextFormat("5.2f");
gStyle->SetOptStat(0);
bool Is_debug=false;
bool kinematiccut= false;
const double DEG=180./3.1410926;   //rad to degree
const double Mp = 938.272081;
const int Nbin=15;
 double bin[Nbin][4]={
 {0.20,0.30,     0.0,14.0},
 {0.30,0.35,     0.0,14.0},
 {0.35,0.40,     0.0, 5.8},
 {0.35,0.40,     5.8,14.0},
 {0.40,0.45,     0.0, 6.4},
 {0.40,0.45,     6.4,14.0},
 {0.45,0.50,     0.0, 7.0},
 {0.45,0.50,     7.0,14.0},
 {0.50,0.55,     0.0, 7.6},
 {0.50,0.55,     7.6,14.0},
 {0.55,0.60,     0.0, 8.2},
 {0.55,0.60,     8.2,14.0},
 {0.60,0.67,     0.0, 8.8},
 {0.60,0.67,     8.8,14.0},
 {0.67,0.80,     0.0,14.0}
 };
double x_cor[15]={0.250,0.325,0.375,0.375,0.425,0.425,0.475,0.475,0.525,0.525,0.575,0.575,0.635,0.635,0.735};
double Q2_cor[15]={4.2,5.0,5.5,6.3,6.0,7.0,6.5,7.8,7.1,8.5,7.6,9.1,8.2,9.8,9.8};
   Double_t x_cor2[12]  = {0.185,0.225,0.265,0.305,0.345,0.385,0.425,0.465,0.505,0.545,0.585,0.625};
   Double_t Q2_cor2[12]  = {2.2,2.8,3.4,4.0,4.6,5.2,5.8,6.4,7.0,7.6,8.2,8.8};
 double cross_LH2_eDIS[12]={ 32.487,16.492,9.356,5.596,3.507,2.229,1.428,0.918,0.588,0.372,0.233,0.142};
 double cross_LH2[12]={ 30.740,15.823,9.024,5.461,3.427,2.221,1.558,1.063,0.718,0.480,0.318,0.206};
 double cross_He3_eDIS[12]={ 90.457,45.066,24.851,14.715,9.103,5.737,3.676,2.329,1.493,0.934,0.569,0.343};
 double cross_He3[12]={ 70.403,35.466,19.854,11.923,7.415,4.680,3.374,2.266,1.518,1.011,0.667,0.434};
 double cross_LD2_eDIS[12]={ 57.704,28.228,15.491,9.072,5.509,3.467,2.176,1.391,0.873,0.548,0.342,0.207};
 double cross_LD2[12]={ 44.877,22.223,12.445,7.348,4.564,2.889,1.988,1.335,0.894,0.596,0.397,0.261};
double ratio_LH2_eDIS[12]={0};
double ratio_He3_eDIS[12]={0};
double ratio_LD2_eDIS[12]={0};
for(int j=0;j<12;j++){
ratio_LH2_eDIS[j]=cross_LH2[j]/cross_LH2_eDIS[j];
ratio_He3_eDIS[j]=cross_He3[j]/cross_He3_eDIS[j];
ratio_LD2_eDIS[j]=cross_LD2[j]/cross_LD2_eDIS[j];
}
TFile *file=new TFile("./PVDIS_LD_interrad_only_486.root");
TFile *file1=new TFile("/w/halla-scifs17exp/solid/tianye/solid_simulation/solid/evgen/eDIS/gen_LD2_5e6_new.root");
TFile *file2=new TFile("./PVDIS_LD2_radoff_847files.root");
TFile *file3=new TFile("./PVDIS_LD2_interrad_Lb_new_485files.root");
TFile *file4=new TFile("./PVDIS_LD2_radall_987files.root");
TFile *file5=new TFile("./PVDIS_LD2_radoff_cetq66_600files.root");
TFile *file6=new TFile("./SIDIS_3He_radall_972files.root");
TFile *file7=new TFile("./SIDIS_3He_norad_976files.root");
TFile *file8=new TFile("/w/halla-scifs17exp/solid/tianye/solid_simulation/solid/evgen/eDIS/gen_He3_5e6_new.root");
TFile *file9=new TFile("/w/halla-scifs17exp/solid/tianye/solid_simulation/solid/evgen/eDIS/gen_LH2_5e6_new.root");
TFile *file10=new TFile("./PVDIS_LH2_norad_125files.root");
TFile *file11=new TFile("./PVDIS_LH2_norad_noscale_300files.root");
TFile *file12=new TFile("./PVDIS_LD2_vertexsmear_Radall_292files.root");
const double filen=486;
const double filen1=987;//418;//492;//484;//515;
const double filen3=485;
const double filenwiser=847;
const double filen_eDIS=1;//5;
const double filen_cetq66=600;
const double filen_SIDIS=972;
const double filen_SIDISno=976;
const double filen_LH=125;
const double filen_LH_noscale=300;
const double filen_LD2_rad=292;


// define histograms, output txt files etc...
        TH1F *p_efnorad[15];
        TH1F *p_efnorad_eDIS[15];
        TH1F *p_efnorad_cetq66[15];
        TH1F *p_efinter[15];
        TH1F *p_efall[15];;
        TH1F *p_efall_rad1[15];;
        TH1F *p_efnorad_rad1[15];;
        TH1F *p_efinter_Lb[15];
        TH2F *Q2_vs_x_norad;
        TH2F *Q2_vs_x_norad_eDIS;
        TH2F *Q2_vs_x_norad_cetq66;
        TH2F *Q2_vs_x_inter;
        TH2F *Q2_vs_x_inter_Lb;
        TH2F *Q2_vs_x_all;
        TH1F *x_norad;
        TH1F *x_norad_cetq66;
        TH1F *x_inter;
        TH1F *x_inter_Lb;
        TH1F *x_all;
        
  char name_norad[50];  
  char name_inter[50];  
  char name_all[50];  
  char name_inter_Lb[50];  
  char name_norad1[50];  
  char name_norad2[50];  
  char name_inter1[50];  
  char name_all1[50];  
  char name_all_rad[50];  
  char name_norad_rad[50];  
  char name_inter_Lb1[50];  
for(int j=0;j<15;j++){
//cout<<"xmin="<<bin[j][0]<<"xmax="<<bin[j][1]<<"Q2min="<<bin[j][2]<<"Q2max="<<bin[j][3]<<endl;
   sprintf(name_norad,"ef_norad_bin%d",j);
   sprintf(name_inter,"ef_inter_bin%d",j);
   sprintf(name_all,"ef_all_bin%d",j);
   sprintf(name_all_rad,"ef_all_rad1_bin%d",j);
   sprintf(name_norad_rad,"ef_norad_rad1_bin%d",j);
   sprintf(name_inter_Lb,"ef_inter_Lb_bin%d",j);
   sprintf(name_norad1,"ef_norad_eDIS_bin%d",j);
   sprintf(name_norad2,"ef_norad_cetq66_bin%d",j);
   p_efnorad[j]=new TH1F(name_norad,name_norad,50,0,10);
   p_efnorad_eDIS[j]=new TH1F(name_norad1,name_norad1,50,0,10);
   p_efnorad_cetq66[j]=new TH1F(name_norad2,name_norad2,50,0,10);
   p_efinter[j]=new TH1F(name_inter,name_inter,50,0,10);
   p_efall[j]=new TH1F(name_all,name_all,50,0,10);
   p_efall_rad1[j]=new TH1F(name_all_rad,name_all_rad,50,0,10);
   p_efnorad_rad1[j]=new TH1F(name_norad_rad,name_norad_rad,50,0,10);
   p_efinter_Lb[j]=new TH1F(name_inter_Lb,name_inter_Lb,50,0,10);
}
   Q2_vs_x_norad=new TH2F("Q2_vs_x_norad","Q2_vs_x_norad",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_xs=new TH2F("Q2_vs_x_norad_xs","Q2_vs_x_norad_xs",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_xs_LH2=new TH2F("Q2_vs_x_norad_xs_LH2","Q2_vs_x_norad_xs_LH2",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_xs_LH2_noscale=new TH2F("Q2_vs_x_norad_xs_LH2_noscale","Q2_vs_x_norad_xs_LH2_noscale",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_W_norad=new TH2F("Q2_vs_W_norad","Q2_vs_W_norad",2,0,6,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_SIDIS=new TH2F("Q2_vs_x_norad_SIDIS","Q2_vs_x_norad_SIDIS",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_LH2=new TH2F("Q2_vs_x_norad_LH2","Q2_vs_x_norad_LH2",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_LH2_noscale=new TH2F("Q2_vs_x_norad_LH2_noscale","Q2_vs_x_norad_LH2_noscale",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_SIDIS=new TH2F("Q2_vs_x_SIDIS","Q2_vs_x_SIDIS",10,0,1,6,0.5,12.5);
   Q2_vs_x_norad_cetq66=new TH2F("Q2_vs_x_norad_cetq66","Q2_vs_x_norad_cetq66",10,0,1,6,0.5,12.5);
   Q2_vs_x_norad_eDIS=new TH2F("Q2_vs_x_norad_eDIS","Q2_vs_x_norad_eDIS",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_eDIS_SIDIS=new TH2F("Q2_vs_x_norad_eDIS_SIDIS","Q2_vs_x_norad_eDIS_SIDIS",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_eDIS_LH2=new TH2F("Q2_vs_x_norad_eDIS_LH2","Q2_vs_x_norad_eDIS_LH2",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_eDIS_xs=new TH2F("Q2_vs_x_norad_eDIS_xs","Q2_vs_x_norad_eDIS_xs",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_eDIS_xs_SIDIS=new TH2F("Q2_vs_x_norad_eDIS_xs_SIDIS","Q2_vs_x_norad_eDIS_xs_SIDIS",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_eDIS_xs_LH2=new TH2F("Q2_vs_x_norad_eDIS_xs_LH2","Q2_vs_x_norad_eDIS_xs_LH2",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_W_norad_eDIS=new TH2F("Q2_vs_W_norad_eDIS","Q2_vs_W_norad_eDIS",2,0,6,6,0.5,12.5);
   Q2_vs_x_inter=new TH2F("Q2_vs_x_inter","Q2_vs_x_internal_rad",10,0,1,6,0.5,12.5);
   Q2_vs_x_inter_Lb=new TH2F("Q2_vs_x_internal_Lb","Q2_vs_x_internal_before_vertex",10,0,1,6,0.5,12.5);
   Q2_vs_x_all=new TH2F("Q2_vs_x_all","Q2_vs_x_radall",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_all_rad1=new TH2F("Q2_vs_x_all_rad1","Q2_vs_x_all_rad1",10,0,1,6,0.5,12.5);
   TH2F *Q2_vs_x_norad_rad1=new TH2F("Q2_vs_x_norad_rad1","Q2_vs_x_norad_rad1",10,0,1,6,0.5,12.5);
   x_norad=new TH1F("x_norad","x_norad",100,0,1);
   x_norad_cetq66=new TH1F("x_norad_cetq66","x_norad_cetq66",100,0,1);
   x_norad_eDIS=new TH1F("x_norad_eDIS","x_norad_eDIS",100,0,1);
   x_inter=new TH1F("x_inter","x_internal_rad",100,0,1);
   x_inter_Lb=new TH1F("x_internal_Lb","x_internal_before_vertex",100,0,1);
   x_all=new TH1F("x_all","x_radall",100,0,1);
   TH2F *Theta_vs_p_norad=new TH2F("Theta_vs_p_norad","Theta_vs_p_norad",10,0,7,10,0,50);
   TH2F *Theta_vs_p_norad_LH2=new TH2F("Theta_vs_p_norad_LH2","Theta_vs_p_norad_LH2",10,0,7,10,0,50);
   TH2F *Theta_vs_p_norad_LH2_noscale=new TH2F("Theta_vs_p_norad_LH2_noscale","Theta_vs_p_norad_LH2_noscale",10,0,7,10,0,50);
   TH2F *Theta_vs_p_cetq66=new TH2F("Theta_vs_p_norad_cetq66","Theta_vs_p_norad_cetq66",10,0,7,10,0,50);
   TH2F *Theta_vs_p_norad_eDIS=new TH2F("Theta_vs_p_norad_eDIS","Theta_vs_p_norad_eDIS",10,0,7,10,0,50);
   TH2F *Theta_vs_p_norad_eDIS_SIDIS=new TH2F("Theta_vs_p_norad_eDIS_SIDIS","Theta_vs_p_norad_eDIS_SIDIS",6,5,35,10,1,11);
   TH2F *Theta_vs_p_norad_eDIS_LH2=new TH2F("Theta_vs_p_norad_eDIS_LH2","Theta_vs_p_norad_eDIS_LH2",10,0,7,10,0,50);
   TH2F *Theta_vs_p_norad_SIDIS=new TH2F("Theta_vs_p_norad_SIDIS","Theta_vs_p_norad_SIDIS",6,5,35,10,1,11);
   TH2F *Theta_vs_p_SIDIS=new TH2F("Theta_vs_p_SIDIS","Theta_vs_p_SIDIS",6,5,35,10,1,11);
   TH2F *Theta_vs_p_inter=new TH2F("Theta_vs_p_inter","Theta_vs_p_internal_rad",10,0,7,10,0,50);
   TH2F *Theta_vs_p_inter_Lb=new TH2F("Theta_vs_p_internal_Lb","Theta_vs_p_internal_before_vertex",10,0,7,10,0,50);
   TH2F *Theta_vs_p_all=new TH2F("Theta_vs_p_all","Theta_vs_p_radall",10,0,7,10,0,50);
   TF1 *f1 = new TF1("f1","(2.0*2.0-0.938*0.938)*x/(1-x)", 0, 1);
   TF1 *f2 = new TF1("f2","(1.23*1.23-0.938*0.938)*x/(1-x)", 0, 1);
   TF1 *f3 = new TF1("f3","(1.48*1.48-0.938*0.938)*x/(1-x)", 0, 1);
   TF1 *f4 = new TF1("f4","(1.66*1.66-0.938*0.938)*x/(1-x)", 0, 1);
   TF1 *f5 = new TF1("f5","(2.3*2.3-0.938*0.938)*x/(1-x)", 0, 1);
   TF1 *f6 = new TF1("f6","(3.0*2.0-0.938*0.938)*x/(1-x)", 0, 1);
	//-------------------------
	//   get trees in the real data file
	//-------------------------
  TLegend* tl = new TLegend (.15, .6, .35, .9);
  tl->AddEntry (f2,"W^{2}=1.23 GeV", "l");
  tl->AddEntry (f3,"W^{2}=1.48 GeV", "l");
  tl->AddEntry (f4,"W^{2}=1.66 GeV", "l");
  tl->AddEntry (f1,"W^{2}=2.0 GeV", "l");
  tl->AddEntry (f6,"W^{2}=3.0 GeV", "l");
//  tl->AddEntry (f5,"W^{2}=2.3 GeV", "l");

f1->SetLineColor(1);
f2->SetLineColor(2);
f3->SetLineColor(6);
f4->SetLineColor(3);
f5->SetLineColor(5);
f6->SetLineColor(7);
f1->SetLineWidth(2);
f2->SetLineWidth(2);
f3->SetLineWidth(2);
f4->SetLineWidth(2);
f5->SetLineWidth(2);
f6->SetLineWidth(2);
	//---header tree
	TTree *tree_header = (TTree*) file->Get("T");
   Double_t        x;
   Double_t        W;
   Double_t        Q2;
   Double_t        rate;
   Double_t        Ep;
   Double_t        vz;
   Double_t        xs;
   Double_t        radxs;
   Double_t        dXSdEdOmega_mubGeVSr;
   Double_t        raddXSdEdOmega_mubGeVSr;
   Double_t        theta;
   Double_t        phi;
   Double_t        nu;

   // List of branches
   TBranch        *b_data;   //!
   tree_header->SetBranchAddress("x", &x, &b_data);
   tree_header->SetBranchAddress("W", &W, &b_data);
   tree_header->SetBranchAddress("Q2", &Q2, &b_data);
   tree_header->SetBranchAddress("rate", &rate, &b_data);
   tree_header->SetBranchAddress("Ep", &Ep, &b_data);
   tree_header->SetBranchAddress("xs", &xs, &b_data);
   tree_header->SetBranchAddress("vz", &vz, &b_data);
   tree_header->SetBranchAddress("radxs", &radxs, &b_data);
   tree_header->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr, &b_data);
   tree_header->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr, &b_data);
   tree_header->SetBranchAddress("theta", &theta, &b_data);
   tree_header->SetBranchAddress("phi", &phi, &b_data);
   tree_header->SetBranchAddress("nu", &nu, &b_data);
	

	TTree *tree_header1 = (TTree*) file4->Get("T");
   Double_t        x1;
   Double_t        W1;
   Double_t        Q21;
   Double_t        rate1;
   Double_t        Ep1;
   Double_t        vz1;
   Double_t        xs1;
   Double_t        radxs1;
   Double_t        dXSdEdOmega_mubGeVSr1;
   Double_t        raddXSdEdOmega_mubGeVSr1;
   Double_t        theta1;
   Double_t        phi1;
   Double_t        nu1;
   tree_header1->SetBranchAddress("x", &x1, &b_data);
   tree_header1->SetBranchAddress("W", &W1, &b_data);
   tree_header1->SetBranchAddress("Q2", &Q21, &b_data);
   tree_header1->SetBranchAddress("rate", &rate1, &b_data);
   tree_header1->SetBranchAddress("Ep", &Ep1, &b_data);
   tree_header1->SetBranchAddress("xs", &xs1, &b_data);
   tree_header1->SetBranchAddress("vz", &vz1, &b_data);
   tree_header1->SetBranchAddress("radxs", &radxs1, &b_data);
   tree_header1->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr1, &b_data);
   tree_header1->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr1, &b_data);
   tree_header1->SetBranchAddress("theta", &theta1, &b_data);
   tree_header1->SetBranchAddress("phi", &phi1, &b_data);
   tree_header1->SetBranchAddress("nu", &nu1, &b_data);


	//---generated tree
	//particle generated with certain momentum at certain vertex
//wiser tree
	TTree *tree_header2 = (TTree*) file2->Get("T");
   Double_t        x2;
   Double_t        W2;
   Double_t        Q22;
   Double_t        rate2;
   Double_t        Ep2;
   Double_t        vz2;
   Double_t        xs2;
   Double_t        radxs2;
   Double_t        dXSdEdOmega_mubGeVSr2;
   Double_t        raddXSdEdOmega_mubGeVSr2;
   Double_t        theta2;
   Double_t        phi2;
   Double_t        nu2;

   // List of branches
   tree_header2->SetBranchAddress("x", &x2, &b_data);
   tree_header2->SetBranchAddress("W", &W2, &b_data);
   tree_header2->SetBranchAddress("Q2", &Q22, &b_data);
   tree_header2->SetBranchAddress("rate", &rate2, &b_data);
   tree_header2->SetBranchAddress("Ep", &Ep2, &b_data);
   tree_header2->SetBranchAddress("xs", &xs2, &b_data);
   tree_header2->SetBranchAddress("vz", &vz2, &b_data);
   tree_header2->SetBranchAddress("radxs", &radxs2, &b_data);
   tree_header2->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr2, &b_data);
   tree_header2->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr2, &b_data);
   tree_header2->SetBranchAddress("theta", &theta2, &b_data);
   tree_header2->SetBranchAddress("phi", &phi2, &b_data);
   tree_header2->SetBranchAddress("nu", &nu2, &b_data);
	
       
	TTree *tree_header3 = (TTree*) file3->Get("T");
   Double_t        x3;
   Double_t        W3;
   Double_t        Q23;
   Double_t        rate3;
   Double_t        Ep3;
   Double_t        vz3;
   Double_t        xs3;
   Double_t        radxs3;
   Double_t        dXSdEdOmega_mubGeVSr3;
   Double_t        raddXSdEdOmega_mubGeVSr3;
   Double_t        theta3;
   Double_t        phi3;
   Double_t        nu3;

   // List of branches
   tree_header3->SetBranchAddress("x", &x3, &b_data);
   tree_header3->SetBranchAddress("W", &W3, &b_data);
   tree_header3->SetBranchAddress("Q2", &Q23, &b_data);
   tree_header3->SetBranchAddress("rate", &rate3, &b_data);
   tree_header3->SetBranchAddress("Ep", &Ep3, &b_data);
   tree_header3->SetBranchAddress("xs", &xs3, &b_data);
   tree_header3->SetBranchAddress("vz", &vz3, &b_data);
   tree_header3->SetBranchAddress("radxs", &radxs3, &b_data);
   tree_header3->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr3, &b_data);
   tree_header3->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr3, &b_data);
   tree_header3->SetBranchAddress("theta", &theta3, &b_data);
   tree_header3->SetBranchAddress("phi", &phi3, &b_data);
   tree_header3->SetBranchAddress("nu", &nu3, &b_data);
	TTree *tree_header4 = (TTree*) file1->Get("T");
   Double_t        weight4;
   TBranch        *b_weight;   //!
   Double_t        rate4;
   Double_t        theta4;
   Double_t        x4;
   Double_t        Ef4;
   Double_t        crs4;
   Double_t        Q24;
   Double_t        W4;
   Double_t        pf4;
   TBranch        *b_rate;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_x;   //!
   TBranch        *b_Ef;   //!
   TBranch        *b_crs;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_pf;   //!
   tree_header4->SetBranchAddress("weight", &weight4, &b_weight);
   tree_header4->SetBranchAddress("rate", &rate4, &b_rate);
   tree_header4->SetBranchAddress("theta", &theta4, &b_theta);
   tree_header4->SetBranchAddress("x", &x4, &b_x);
   tree_header4->SetBranchAddress("Ef", &Ef4, &b_Ef);
   tree_header4->SetBranchAddress("crs", &crs4, &b_crs);
   tree_header4->SetBranchAddress("Q2", &Q24, &b_Q2);
   tree_header4->SetBranchAddress("W", &W4, &b_W);
   tree_header4->SetBranchAddress("pf", &pf4, &b_pf);
	TTree *tree_header5 = (TTree*) file5->Get("T");
   Double_t        x5;
   Double_t        W5;
   Double_t        Q25;
   Double_t        rate5;
   Double_t        Ep5;
   Double_t        vz5;
   Double_t        xs5;
   Double_t        radxs5;
   Double_t        dXSdEdOmega_mubGeVSr5;
   Double_t        raddXSdEdOmega_mubGeVSr5;
   Double_t        theta5;
   Double_t        phi5;
   Double_t        nu5;

   // List of branches
   tree_header5->SetBranchAddress("x", &x5, &b_data);
   tree_header5->SetBranchAddress("W", &W5, &b_data);
   tree_header5->SetBranchAddress("Q2", &Q25, &b_data);
   tree_header5->SetBranchAddress("rate", &rate5, &b_data);
   tree_header5->SetBranchAddress("Ep", &Ep5, &b_data);
   tree_header5->SetBranchAddress("xs", &xs5, &b_data);
   tree_header5->SetBranchAddress("vz", &vz5, &b_data);
   tree_header5->SetBranchAddress("radxs", &radxs5, &b_data);
   tree_header5->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr5, &b_data);
   tree_header5->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr5, &b_data);
   tree_header5->SetBranchAddress("theta", &theta5, &b_data);
   tree_header5->SetBranchAddress("phi", &phi5, &b_data);
   tree_header5->SetBranchAddress("nu", &nu5, &b_data);
	TTree *tree_header6 = (TTree*) file6->Get("T");
   Double_t        x6;
   Double_t        W6;
   Double_t        Q26;
   Double_t        rate6;
   Double_t        Ep6;
   Double_t        vz6;
   Double_t        xs6;
   Double_t        radxs6;
   Double_t        dXSdEdOmega_mubGeVSr6;
   Double_t        raddXSdEdOmega_mubGeVSr6;
   Double_t        theta6;
   Double_t        phi6;
   Double_t        nu6;

   // List of branches
   TBranch        *b_data;   //!
   tree_header6->SetBranchAddress("x", &x6, &b_data);
   tree_header6->SetBranchAddress("W", &W6, &b_data);
   tree_header6->SetBranchAddress("Q2", &Q26, &b_data);
   tree_header6->SetBranchAddress("rate", &rate6, &b_data);
   tree_header6->SetBranchAddress("Ep", &Ep6, &b_data);
   tree_header6->SetBranchAddress("xs", &xs6, &b_data);
   tree_header6->SetBranchAddress("vz", &vz6, &b_data);
   tree_header6->SetBranchAddress("radxs", &radxs6, &b_data);
   tree_header6->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr6, &b_data);
   tree_header6->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr6, &b_data);
   tree_header6->SetBranchAddress("theta", &theta6, &b_data);
   tree_header6->SetBranchAddress("phi", &phi6, &b_data);
   tree_header6->SetBranchAddress("nu", &nu6, &b_data);
	TTree *tree_header7 = (TTree*) file7->Get("T");
   Double_t        x7;
   Double_t        W7;
   Double_t        Q27;
   Double_t        rate7;
   Double_t        Ep7;
   Double_t        vz7;
   Double_t        xs7;
   Double_t        radxs7;
   Double_t        dXSdEdOmega_mubGeVSr7;
   Double_t        raddXSdEdOmega_mubGeVSr7;
   Double_t        theta7;
   Double_t        phi7;
   Double_t        nu7;

   // List of branches
   TBranch        *b_data;   //!
   tree_header7->SetBranchAddress("x", &x7, &b_data);
   tree_header7->SetBranchAddress("W", &W7, &b_data);
   tree_header7->SetBranchAddress("Q2", &Q27, &b_data);
   tree_header7->SetBranchAddress("rate", &rate7, &b_data);
   tree_header7->SetBranchAddress("Ep", &Ep7, &b_data);
   tree_header7->SetBranchAddress("xs", &xs7, &b_data);
   tree_header7->SetBranchAddress("vz", &vz7, &b_data);
   tree_header7->SetBranchAddress("radxs", &radxs7, &b_data);
   tree_header7->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr7, &b_data);
   tree_header7->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr7, &b_data);
   tree_header7->SetBranchAddress("theta", &theta7, &b_data);
   tree_header7->SetBranchAddress("phi", &phi7, &b_data);
   tree_header7->SetBranchAddress("nu", &nu7, &b_data);
	
	TTree *tree_header8 = (TTree*) file8->Get("T");
   Double_t        weight8;
   TBranch        *b_weight;   //!
   Double_t        rate8;
   Double_t        theta8;
   Double_t        x8;
   Double_t        Ef8;
   Double_t        crs8;
   Double_t        Q28;
   Double_t        W8;
   Double_t        pf8;
   TBranch        *b_rate;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_x;   //!
   TBranch        *b_Ef;   //!
   TBranch        *b_crs;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_pf;   //!
   tree_header8->SetBranchAddress("weight", &weight8, &b_weight);
   tree_header8->SetBranchAddress("rate", &rate8, &b_rate);
   tree_header8->SetBranchAddress("theta", &theta8, &b_theta);
   tree_header8->SetBranchAddress("x", &x8, &b_x);
   tree_header8->SetBranchAddress("Ef", &Ef8, &b_Ef);
   tree_header8->SetBranchAddress("crs", &crs8, &b_crs);
   tree_header8->SetBranchAddress("Q2", &Q28, &b_Q2);
   tree_header8->SetBranchAddress("W", &W8, &b_W);
   tree_header8->SetBranchAddress("pf", &pf8, &b_pf);
	TTree *tree_header9 = (TTree*) file9->Get("T");
   Double_t        weight9;
   TBranch        *b_weight;   //!
   Double_t        rate9;
   Double_t        theta9;
   Double_t        x9;
   Double_t        Ef9;
   Double_t        crs9;
   Double_t        Q29;
   Double_t        W9;
   Double_t        pf9;
   TBranch        *b_rate;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_x;   //!
   TBranch        *b_Ef;   //!
   TBranch        *b_crs;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_pf;   //!
   tree_header9->SetBranchAddress("weight", &weight9, &b_weight);
   tree_header9->SetBranchAddress("rate", &rate9, &b_rate);
   tree_header9->SetBranchAddress("theta", &theta9, &b_theta);
   tree_header9->SetBranchAddress("x", &x9, &b_x);
   tree_header9->SetBranchAddress("Ef", &Ef9, &b_Ef);
   tree_header9->SetBranchAddress("crs", &crs9, &b_crs);
   tree_header9->SetBranchAddress("Q2", &Q29, &b_Q2);
   tree_header9->SetBranchAddress("W", &W9, &b_W);
   tree_header9->SetBranchAddress("pf", &pf9, &b_pf);
	TTree *tree_header10 = (TTree*) file10->Get("T");
   Double_t        x10;
   Double_t        W10;
   Double_t        Q210;
   Double_t        rate10;
   Double_t        Ep10;
   Double_t        vz10;
   Double_t        xs10;
   Double_t        radxs10;
   Double_t        dXSdEdOmega_mubGeVSr10;
   Double_t        raddXSdEdOmega_mubGeVSr10;
   Double_t        theta10;
   Double_t        phi10;
   Double_t        nu10;

   // List of branches
   TBranch        *b_data;   //!
   tree_header10->SetBranchAddress("x", &x10, &b_data);
   tree_header10->SetBranchAddress("W", &W10, &b_data);
   tree_header10->SetBranchAddress("Q2", &Q210, &b_data);
   tree_header10->SetBranchAddress("rate", &rate10, &b_data);
   tree_header10->SetBranchAddress("Ep", &Ep10, &b_data);
   tree_header10->SetBranchAddress("xs", &xs10, &b_data);
   tree_header10->SetBranchAddress("vz", &vz10, &b_data);
   tree_header10->SetBranchAddress("radxs", &radxs10, &b_data);
   tree_header10->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr10, &b_data);
   tree_header10->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr10, &b_data);
   tree_header10->SetBranchAddress("theta", &theta10, &b_data);
   tree_header10->SetBranchAddress("phi", &phi10, &b_data);
   tree_header10->SetBranchAddress("nu", &nu10, &b_data);
	
	TTree *tree_header11 = (TTree*) file11->Get("T");
   Double_t        x11;
   Double_t        W11;
   Double_t        Q211;
   Double_t        rate11;
   Double_t        Ep11;
   Double_t        vz11;
   Double_t        xs11;
   Double_t        radxs11;
   Double_t        dXSdEdOmega_mubGeVSr11;
   Double_t        raddXSdEdOmega_mubGeVSr11;
   Double_t        theta11;
   Double_t        phi11;
   Double_t        nu11;

   // List of branches
   TBranch        *b_data;   //!
   tree_header11->SetBranchAddress("x", &x11, &b_data);
   tree_header11->SetBranchAddress("W", &W11, &b_data);
   tree_header11->SetBranchAddress("Q2", &Q211, &b_data);
   tree_header11->SetBranchAddress("rate", &rate11, &b_data);
   tree_header11->SetBranchAddress("Ep", &Ep11, &b_data);
   tree_header11->SetBranchAddress("xs", &xs11, &b_data);
   tree_header11->SetBranchAddress("vz", &vz11, &b_data);
   tree_header11->SetBranchAddress("radxs", &radxs11, &b_data);
   tree_header11->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr11, &b_data);
   tree_header11->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr11, &b_data);
   tree_header11->SetBranchAddress("theta", &theta11, &b_data);
   tree_header11->SetBranchAddress("phi", &phi11, &b_data);
   tree_header11->SetBranchAddress("nu", &nu11, &b_data);
   TTree *tree_header12 = (TTree*) file12->Get("T");
   Double_t        x12;
   Double_t        W12;
   Double_t        Q212;
   Double_t        rate12;
   Double_t        raterad12;
   Double_t        Ep12;
   Double_t        vz12;
   Double_t        xs12;
   Double_t        radxs12;
   Double_t        dXSdEdOmega_mubGeVSr12;
   Double_t        raddXSdEdOmega_mubGeVSr12;
   Double_t        theta12;
   Double_t        phi12;
   Double_t        nu12;
   tree_header12->SetBranchAddress("x", &x12, &b_data);
   tree_header12->SetBranchAddress("W", &W12, &b_data);
   tree_header12->SetBranchAddress("Q2", &Q212, &b_data);
   tree_header12->SetBranchAddress("rate", &rate12, &b_data);
   tree_header12->SetBranchAddress("raterad", &raterad12, &b_data);
   tree_header12->SetBranchAddress("Ep", &Ep12, &b_data);
   tree_header12->SetBranchAddress("xs", &xs12, &b_data);
   tree_header12->SetBranchAddress("vz", &vz12, &b_data);
   tree_header12->SetBranchAddress("radxs", &radxs12, &b_data);
   tree_header12->SetBranchAddress("dXSdEdOmega_mubGeVSr", &dXSdEdOmega_mubGeVSr12, &b_data);
   tree_header12->SetBranchAddress("raddXSdEdOmega_mubGeVSr", &raddXSdEdOmega_mubGeVSr12, &b_data);
   tree_header12->SetBranchAddress("theta", &theta12, &b_data);
   tree_header12->SetBranchAddress("phi", &phi12, &b_data);
   tree_header12->SetBranchAddress("nu", &nu12, &b_data);

	long int N_eventsw = (long int)tree_header2->GetEntries();
	long int N_events3 = (long int)tree_header3->GetEntries();
	long int N_events = (long int)tree_header->GetEntries();
	long int N_events1 = (long int)tree_header1->GetEntries();
	long int N_events4 = (long int)tree_header4->GetEntries();
	long int N_events5 = (long int)tree_header5->GetEntries();
	long int N_events6 = (long int)tree_header6->GetEntries();
	long int N_events7 = (long int)tree_header7->GetEntries();
	long int N_events8 = (long int)tree_header8->GetEntries();
	long int N_events9 = (long int)tree_header9->GetEntries();
	long int N_events10 = (long int)tree_header10->GetEntries();
	long int N_events11 = (long int)tree_header11->GetEntries();
	long int N_events12 = (long int)tree_header12->GetEntries();
        cout<<"N_inclusive="<<N_eventsw<<"N_eDIS="<<N_events<<endl;
	for(long int i=0;i<N_events;i++){	  
		tree_header->GetEntry(i);
		double ratef=rate/filen;
              //  cout<<"rate="<<rate<<endl;
                double theta_eDIS=theta;
              if(rate>0 /*&& W>2*/){
              Q2_vs_x_inter->Fill(x, Q2, ratef); 
              Theta_vs_p_inter->Fill(Ep, theta_eDIS, ratef); 
              x_inter->Fill(x,  ratef); 
              for(int ibin=0;ibin<Nbin;ibin++){
                if(W>2 && x>bin[ibin][0] && x<bin[ibin][1] && Q2>bin[ibin][2] && Q2<bin[ibin][3]){
                    p_efinter[ibin]->Fill(Ep,ratef);
                }
               }
              }      
         }//end event loop
	for(long int i=0;i<N_events1;i++){	  
		tree_header1->GetEntry(i);
		double ratef1=rate1/filen1;
                double theta_eDIS1=theta1; 
              if(rate1>0 /*&& W1>2*/){
              Q2_vs_x_all->Fill(x1, Q21, ratef1); 
              Theta_vs_p_all->Fill(Ep1, theta_eDIS1, ratef1); 
              x_all->Fill(x1, ratef1); 
              for(int ibin=0;ibin<Nbin;ibin++){
                if(W1>2 && x1>bin[ibin][0] && x1<bin[ibin][1] && Q21>bin[ibin][2] && Q21<bin[ibin][3]){
                    p_efall[ibin]->Fill(Ep1,ratef1);
                }
              }      
           }
         }
	for(long int i=0;i<N_eventsw;i++){	  
		//---header tree
		//---
		tree_header2->GetEntry(i);
		double ratef2=rate2/filenwiser;
                double  p_gen=Ep2;
                double  theta_gen=theta2;
                double  cross2=dXSdEdOmega_mubGeVSr2*1000;
//                cout<<"cross2="<<cross2<<endl;
              if(rate2>0 /*&& W2>2*/){
              Theta_vs_p_norad->Fill(Ep2, theta_gen, ratef2); 
              Q2_vs_x_norad->Fill(x2, Q22, ratef2); 
              Q2_vs_x_norad_xs->Fill(x2, Q22,cross2); 
              Q2_vs_W_norad->Fill(W2, Q22, ratef2); 
              x_norad->Fill(x2, ratef2); 
              for(int ibin=0;ibin<Nbin;ibin++){
                if(W2>2 && x2>bin[ibin][0] && x2<bin[ibin][1] && Q22>bin[ibin][2] && Q22<bin[ibin][3]){
                    p_efnorad[ibin]->Fill(Ep2,ratef2);
                }
              }  
            }    
         }//end event loop
	for(long int i=0;i<N_events3;i++){	  
		//---header tree
		//---
		tree_header3->GetEntry(i);
		double ratef3=rate3/filen3;
                double  p_gen3=Ep3;
                double  theta_gen3=theta3;
              if(rate3>0  /*&& W3>2*/){
              Theta_vs_p_inter_Lb->Fill(Ep3, theta_gen3, ratef3); 
              Q2_vs_x_inter_Lb->Fill(x3, Q23, ratef3); 
              x_inter_Lb->Fill(x3, ratef3); 
              for(int ibin=0;ibin<Nbin;ibin++){
                if(W3>2 && x3>bin[ibin][0] && x3<bin[ibin][1] && Q23>bin[ibin][2] && Q23<bin[ibin][3]){
                    p_efinter_Lb[ibin]->Fill(Ep3,ratef3);
                }
              }  
           }    
         }//end event loop
	for(long int i=0;i<N_events4;i++){	  
		tree_header4->GetEntry(i);
		double rate_eDIS=rate4/filen_eDIS;
                double theta_eDIS4=theta4*DEG; 
                double cross4 = crs4/filen_eDIS;
                //cout<<"cross4="<<cross4<<endl;
              if(rate4>0  /*&& W4>2*/){
              Theta_vs_p_norad_eDIS->Fill(pf4, theta_eDIS4, rate_eDIS); 
              Q2_vs_x_norad_eDIS->Fill(x4, Q24, rate_eDIS); 
              Q2_vs_x_norad_eDIS_xs->Fill(x4, Q24,cross4); 
              Q2_vs_W_norad_eDIS->Fill(W4, Q24, rate_eDIS); 
              x_norad_eDIS->Fill(x4, rate_eDIS); 
              for(int ibin=0;ibin<Nbin;ibin++){
                if(W4>2 && x4>bin[ibin][0] && x4<bin[ibin][1] && Q24>bin[ibin][2] && Q24<bin[ibin][3]){
                    p_efnorad_eDIS[ibin]->Fill(pf4,rate_eDIS);
                }
              }  
           }    
         }//end event loop
	for(long int i=0;i<N_events8;i++){	  
		tree_header8->GetEntry(i);
		double rate_eDIS8=rate8/filen_eDIS;
                double theta_eDIS8=theta8*DEG; 
                double cross8 = crs8/filen_eDIS;
                //cout<<"cross4="<<cross4<<endl;
              if(rate8>0  ){
              Theta_vs_p_norad_eDIS_SIDIS->Fill( theta_eDIS8,pf8 ,rate_eDIS8); 
              Q2_vs_x_norad_eDIS_SIDIS->Fill(x8, Q28, rate_eDIS8); 
              Q2_vs_x_norad_eDIS_xs_SIDIS->Fill(x8, Q28,cross8); 
           }    
         }//end event loop
	for(long int i=0;i<N_events9;i++){	  
		tree_header9->GetEntry(i);
		double rate_eDIS9=rate9/filen_eDIS;
                double theta_eDIS9=theta9*DEG; 
                double cross9 = crs9/filen_eDIS;
                //cout<<"cross4="<<cross4<<endl;
              if(rate9>0  ){
              Theta_vs_p_norad_eDIS_LH2->Fill(pf9, theta_eDIS9, rate_eDIS9); 
              Q2_vs_x_norad_eDIS_LH2->Fill(x9, Q29, rate_eDIS9); 
              Q2_vs_x_norad_eDIS_xs_LH2->Fill(x9, Q29,cross9); 
           }    
         }//end event loop
	for(long int i=0;i<N_events5;i++){	  
		tree_header5->GetEntry(i);
		double rate66=rate5/filen_cetq66;
                double theta_eDIS5=theta5; 
              if(rate5>0  /*&& W5>2*/){
              Theta_vs_p_norad_cetq66->Fill(Ep5, theta_eDIS5, rate66); 
              Q2_vs_x_norad_cetq66->Fill(x5, Q25, rate66); 
              x_norad_cetq66->Fill(x5, rate66); 
              for(int ibin=0;ibin<Nbin;ibin++){
                if(W5>2 && x5>bin[ibin][0] && x5<bin[ibin][1] && Q25>bin[ibin][2] && Q25<bin[ibin][3]){
                    p_efnorad_cetq66[ibin]->Fill(Ep5,rate66);
                }
              }  
           }
        }    
	for(long int i=0;i<N_events6;i++){	  
		tree_header6->GetEntry(i);
		double ratef6=rate6/filen_SIDIS;
                double theta_eDIS6=theta6; 
              if(ratef6>0 ){
              Theta_vs_p_SIDIS->Fill(theta_eDIS6,Ep6, ratef6); 
              Q2_vs_x_SIDIS->Fill(x6, Q26, ratef6); 
           }
        }    
	for(long int i=0;i<N_events7;i++){	  
		tree_header7->GetEntry(i);
		double ratef7=rate7/filen_SIDISno;
                double theta_eDIS7=theta7; 
              if(ratef7>0 ){
              Theta_vs_p_norad_SIDIS->Fill(theta_eDIS7,Ep7, ratef7); 
              Q2_vs_x_norad_SIDIS->Fill(x7, Q27, ratef7); 
              }  
         }//end event loop
	for(long int i=0;i<N_events10;i++){	  
		tree_header10->GetEntry(i);
		double ratef10=rate10/filen_LH;
                double theta_eDIS10=theta10; 
                double  cross10=dXSdEdOmega_mubGeVSr10*1000;
              if(ratef10>0  ){
              Theta_vs_p_norad_LH2->Fill(Ep10, theta_eDIS10, ratef10); 
              Q2_vs_x_norad_LH2->Fill(x10, Q210, ratef10); 
              Q2_vs_x_norad_xs_LH2->Fill(x10, Q210, cross10); 
              }  
         }//end event loop
	for(long int i=0;i<N_events11;i++){	  
		tree_header11->GetEntry(i);
		double ratef11=rate11/filen_LH_noscale;
                double theta_eDIS11=theta11; 
                double  cross11=dXSdEdOmega_mubGeVSr11*1000;
              if(ratef11>0  ){
              Theta_vs_p_norad_LH2_noscale->Fill(Ep11, theta_eDIS11, ratef11); 
              Q2_vs_x_norad_LH2_noscale->Fill(x11, Q211, ratef11); 
              Q2_vs_x_norad_xs_LH2_noscale->Fill(x11, Q211, cross11); 
              }  
         }//end event loop
	for(long int i=0;i<N_events12;i++){	  
		tree_header12->GetEntry(i);
		double ratef12=rate12/filen_LD2_rad;
		double raterad12=raterad12/filen_LD2_rad;
                double theta_eDIS12=theta12; 
              if(ratef12>0  ){
              Q2_vs_x_all_rad1->Fill(x12, Q212, raterad12); 
              Q2_vs_x_norad_rad1->Fill(x12, Q212, ratef12); 
              for(int ibin=0;ibin<Nbin;ibin++){
                if(W12>2 && x12>bin[ibin][0] && x12<bin[ibin][1] && Q212>bin[ibin][2] && Q212<bin[ibin][3]){
                    p_efall_rad1[ibin]->Fill(Ep12,raterad12);
                    p_efnorad_rad1[ibin]->Fill(Ep12,rate12);
                }
              }      
              }  
         }//end event loop
         //wiser kl0
        double ratio_rad_onoff[15];
        double ratio_interrad_off[15];
        double ratio_interrad_Lb_off[15];
        double ratio_eDIS_eAll[15];
        double ratio_rad_onoff_rad[15];
        double totalN[15];
        double totalN_norad[15];
  for(int k=0;k<Nbin;k++){
        totalN[k]=p_efall[k]->Integral()/1000.0;
        totalN_norad[k]=p_efnorad[k]->Integral()/1000.0;
        ratio_rad_onoff[k]=p_efall[k]->Integral()/p_efnorad[k]->Integral(); 
        ratio_rad_onoff_rad[k]=p_efall_rad1[k]->Integral()/p_efnorad_rad1[k]->Integral(); 
        ratio_eDIS_eAll[k]=p_efnorad[k]->Integral()/p_efnorad_eDIS[k]->Integral(); 
        ratio_interrad_off[k]=p_efinter[k]->Integral()/p_efnorad[k]->Integral(); 
        ratio_interrad_Lb_off[k]=p_efinter_Lb[k]->Integral()/p_efnorad[k]->Integral(); 
   }  
       TH2F *Q2_vs_x_all_halftarget = (TH2F*) Q2_vs_x_all->Clone();  
       Q2_vs_x_all_halftarget->SetName("Q2_vs_x_all_halftarget");
        TCanvas* c111;
        c111 = new TCanvas ("c111", " Q2_vs_x_rad", 800, 800);
        Q2_vs_x_all_rad1->Divide(Q2_vs_x_norad_rad1);       
        Q2_vs_x_all_rad1->GetXaxis()->SetTitle("x");       
        Q2_vs_x_all_rad1->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_all_rad1->SetTitle("Rate_{rad_vertex}/Rate_{norad}");       
        Q2_vs_x_all_rad1->Draw("COLZ text");
        Q2_vs_x_all_rad1->SetMarkerSize(2);
        Q2_vs_x_all_rad1->SetMarkerColor(1);
for(int k = 0; k < Nbin; k++){
TMarker marker12;
marker12.SetMarkerStyle(20);
marker12.SetMarkerColor(kRed);
marker12.DrawMarker(x_cor[k],Q2_cor[k]);
TText *label12 = new TText(x_cor[k],Q2_cor[k],Form("%.02f",ratio_rad_onoff_rad[k]));
label12->SetTextColor(2);
label12->SetTextSize(0.04);
label12->Draw();
}
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
        TCanvas* c11;
        c11 = new TCanvas ("c11", " Q2_vs_x", 800, 800);
        Q2_vs_x_all->Divide(Q2_vs_x_norad);       
        Q2_vs_x_all->GetXaxis()->SetTitle("x");       
        Q2_vs_x_all->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_all->SetTitle("Rate_{rad}/Rate_{norad}");       
        Q2_vs_x_all->Draw("COLZ text");
        Q2_vs_x_all->SetMarkerSize(2);
        Q2_vs_x_all->SetMarkerColor(1);
for(int k = 0; k < Nbin; k++){
TMarker marker;
marker.SetMarkerStyle(20);
marker.SetMarkerColor(kRed);
marker.DrawMarker(x_cor[k],Q2_cor[k]);
TText *label = new TText(x_cor[k],Q2_cor[k],Form("%.02f",ratio_rad_onoff[k]));
label->SetTextColor(2);
label->SetTextSize(0.04);
label->Draw();
}
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
        TCanvas* c13;
        c13 = new TCanvas ("c13", " Q2_vs_x", 800, 800);
        Q2_vs_x_inter_Lb->Divide(Q2_vs_x_norad);       
        Q2_vs_x_inter_Lb->GetXaxis()->SetTitle("x");       
        Q2_vs_x_inter_Lb->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_inter_Lb->SetTitle("Rate_{internal+Lb}/Rate_{norad}");       
        Q2_vs_x_inter_Lb->Draw("COLZ text");
        Q2_vs_x_inter_Lb->SetMarkerSize(2);
        Q2_vs_x_inter_Lb->SetMarkerColor(1);
for(int k = 0; k < Nbin; k++){
TMarker marker1;
marker1.SetMarkerStyle(20);
marker1.SetMarkerColor(kRed);
marker1.DrawMarker(x_cor[k],Q2_cor[k]);
TText *label1 = new TText(x_cor[k],Q2_cor[k],Form("%.02f",ratio_interrad_off[k]));
label1->SetTextColor(2);
label1->SetTextSize(0.04);
label1->Draw();
}
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
        TCanvas* c12;
        c12 = new TCanvas ("c12", " x", 800, 800);
        x_all->SetLineColor(1);
        x_all->Draw();
        x_inter->SetLineColor(2);
        x_inter->Draw("same");
        x_inter_Lb->SetLineColor(4);
        x_inter_Lb->Draw("same");
        x_norad->SetLineColor(6);
        x_norad->Draw("same");
        x_norad_eDIS->SetLineColor(7);
        x_norad_eDIS->Draw("same");

        TCanvas* c2;
        c2 = new TCanvas ("c2", " theta_vs_p", 800, 800);
        Theta_vs_p_all->Divide(Theta_vs_p_norad);       
        Theta_vs_p_all->Draw("COLZ text");
        Theta_vs_p_all->SetMarkerSize(2.5);
        Theta_vs_p_all->SetMarkerColor(0);
        TCanvas* c3;
        c3 = new TCanvas ("c3", " Q2_vs_x_cetq66", 800, 800);
        Q2_vs_x_norad_cetq66->Divide(Q2_vs_x_norad);       
        Q2_vs_x_norad_cetq66->GetXaxis()->SetTitle("x");       
        Q2_vs_x_norad_cetq66->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_norad_cetq66->SetTitle("Rate_{cteq66}/Rate_{CT14nlo}");       
        Q2_vs_x_norad_cetq66->Draw("COLZ text");
        Q2_vs_x_norad_cetq66->SetMarkerSize(2.5);
        Q2_vs_x_norad_cetq66->SetMarkerColor(0);
        TCanvas* c1;
        c1 = new TCanvas ("c1", " Q2_vs_x", 800, 800);
        Q2_vs_x_norad->Divide(Q2_vs_x_norad_eDIS);       
        Q2_vs_x_norad->GetXaxis()->SetTitle("x");       
        Q2_vs_x_norad->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_norad->SetTitle("Rate_{eAll}/Rate_{eDIS}");       
        Q2_vs_x_norad->Draw("COLZ text");
        Q2_vs_x_norad->SetMarkerSize(2.5);
        Q2_vs_x_norad->SetMarkerColor(0);
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
f6->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
for(int k = 0; k < Nbin; k++){
TMarker marker2;
marker2.SetMarkerStyle(20);
marker2.SetMarkerColor(kRed);
marker2.DrawMarker(x_cor[k],Q2_cor[k]);
TText *label2 = new TText(x_cor[k],Q2_cor[k],Form("%.02f",ratio_eDIS_eAll[k]));
label2->SetTextColor(2);
label2->SetTextSize(0.04);
label2->Draw();
}
        TCanvas* c14;
        c14 = new TCanvas ("c14", " Q2_vs_W", 800, 800);
        Q2_vs_W_norad->Divide(Q2_vs_W_norad_eDIS);       
        Q2_vs_W_norad->GetXaxis()->SetTitle("W [GeV]");       
        Q2_vs_W_norad->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_W_norad->SetTitle("Rate_{eAll}/Rate_{eDIS}");       
        Q2_vs_W_norad->Draw("COLZ text");
        Q2_vs_W_norad->SetMarkerSize(2.5);
        Q2_vs_W_norad->SetMarkerColor(0);
        TCanvas* c15;
        c15 = new TCanvas ("c15", " Q2_vs_x_xs", 800, 800);
        Q2_vs_x_norad_xs->Divide(Q2_vs_x_norad_eDIS_xs);       
        Q2_vs_x_norad_xs->GetXaxis()->SetTitle("W [GeV]");       
        Q2_vs_x_norad_xs->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_norad_xs->SetTitle("#sigma_{eAll}/#sigma_{eDIS}");       
        Q2_vs_x_norad_xs->Draw("COLZ text");
        Q2_vs_x_norad_xs->SetMarkerSize(2.5);
        Q2_vs_x_norad_xs->SetMarkerColor(0);
        TCanvas* c4;
        c4 = new TCanvas ("c4", " Q2_vs_x", 800, 800);
        Q2_vs_x_SIDIS->Divide(Q2_vs_x_norad_SIDIS);       
        Q2_vs_x_SIDIS->GetXaxis()->SetTitle("x");       
        Q2_vs_x_SIDIS->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_SIDIS->SetTitle("Rate_{rad_SIDIS}/Rate_{norad_SIDIS}");       
        Q2_vs_x_SIDIS->Draw("COLZ text");
        Q2_vs_x_SIDIS->SetMarkerSize(2);
        Q2_vs_x_SIDIS->SetMarkerColor(1);
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
f6->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
        TCanvas* c5;
        c5 = new TCanvas ("c5", " Q2_vs_x_SIDIS", 800, 800);
        Q2_vs_x_norad_SIDIS->Divide(Q2_vs_x_norad_eDIS_SIDIS);       
        Q2_vs_x_norad_SIDIS->GetXaxis()->SetTitle("x");       
        Q2_vs_x_norad_SIDIS->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_norad_SIDIS->SetTitle("Rate^{SIDIS}_{eAll}/Rate^{SIDIS}_{eDIS}");       
        Q2_vs_x_norad_SIDIS->Draw("COLZ text");
        Q2_vs_x_norad_SIDIS->SetMarkerSize(2.5);
        Q2_vs_x_norad_SIDIS->SetMarkerColor(0);
for(int k = 0; k < 12; k++){
TMarker marker10;
marker10.SetMarkerStyle(20);
marker10.SetMarkerColor(kRed);
marker10.DrawMarker(x_cor2[k],Q2_cor2[k]);
TText *label10 = new TText(x_cor2[k],Q2_cor2[k],Form("%.02f",ratio_He3_eDIS[k]));
label10->SetTextColor(2);
label10->SetTextSize(0.04);
label10->Draw();
}
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
f6->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
        TCanvas* c6;
        c6 = new TCanvas ("c6", " Q2_vs_x_PVDIS", 800, 800);
        Q2_vs_x_norad_LH2->Divide(Q2_vs_x_norad_eDIS_LH2);       
        Q2_vs_x_norad_LH2->GetXaxis()->SetTitle("x");       
        Q2_vs_x_norad_LH2->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_norad_LH2->SetTitle("Rate^{LH2}_{eAll}/Rate^{LH2}_{eDIS}");       
        Q2_vs_x_norad_LH2->Draw("COLZ text");
        Q2_vs_x_norad_LH2->SetMarkerSize(2.5);
        Q2_vs_x_norad_LH2->SetMarkerColor(0);
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
f6->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
        TCanvas* c7;
        c7 = new TCanvas ("c7", " Q2_vs_x_PVDIS noscale", 800, 800);
        Q2_vs_x_norad_LH2_noscale->Divide(Q2_vs_x_norad_eDIS_LH2);       
        Q2_vs_x_norad_LH2_noscale->GetXaxis()->SetTitle("x");       
        Q2_vs_x_norad_LH2_noscale->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_norad_LH2_noscale->SetTitle("Rate^{LH2 noscale}_{eAll}/Rate^{LH2}_{eDIS}");       
        Q2_vs_x_norad_LH2_noscale->Draw("COLZ text");
        Q2_vs_x_norad_LH2_noscale->SetMarkerSize(2.5);
        Q2_vs_x_norad_LH2_noscale->SetMarkerColor(0);
for(int k = 0; k < 12; k++){
TMarker marker0;
marker0.SetMarkerStyle(20);
marker0.SetMarkerColor(kRed);
marker0.DrawMarker(x_cor2[k],Q2_cor2[k]);
TText *label0 = new TText(x_cor2[k],Q2_cor2[k],Form("%.02f",ratio_LH2_eDIS[k]));
label0->SetTextColor(2);
label0->SetTextSize(0.04);
label0->Draw();
}
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
f6->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
        TCanvas* c8;
        c8 = new TCanvas ("c8", " Q2_vs_x_PVDIS_LH2", 800, 800);
        Q2_vs_x_norad_xs_LH2_noscale->Divide(Q2_vs_x_norad_eDIS_xs_LH2);       
        Q2_vs_x_norad_xs_LH2_noscale->GetXaxis()->SetTitle("x");       
        Q2_vs_x_norad_xs_LH2_noscale->GetYaxis()->SetTitle("Q2 [GeV^{2}]");       
        Q2_vs_x_norad_xs_LH2_noscale->SetTitle("#sigma^{LH2}_{eAll}/#sigma^{LH2}_{eDIS}");       
        Q2_vs_x_norad_xs_LH2_noscale->Draw("COLZ text");
        Q2_vs_x_norad_xs_LH2_noscale->SetMarkerSize(2.5);
        Q2_vs_x_norad_xs_LH2_noscale->SetMarkerColor(0);
f1->Draw("same");
f2->Draw("same");
f3->Draw("same");
f4->Draw("same");
f6->Draw("same");
  tl->SetTextSize (0.015);
  tl->Draw();
        TCanvas* c9;
        c9 = new TCanvas ("c9", " theta_vs_p_SIDIS", 800, 800);
        Theta_vs_p_norad_SIDIS->Divide(Theta_vs_p_norad_eDIS_SIDIS);       
        Theta_vs_p_norad_SIDIS->GetXaxis()->SetTitle("#theta");       
        Theta_vs_p_norad_SIDIS->GetYaxis()->SetTitle("p [GeV]");       
        Theta_vs_p_norad_SIDIS->SetTitle("Rate^{SIDIS}_{eAll}/Rate^{SIDIS}_{eDIS}");       
        Theta_vs_p_norad_SIDIS->Draw("COLZ text");
        Theta_vs_p_norad_SIDIS->SetMarkerSize(2.5);
        Theta_vs_p_norad_SIDIS->SetMarkerColor(0);
} 

