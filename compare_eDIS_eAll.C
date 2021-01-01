// test script (?) 

#include <cstdlib>
#include <iostream> 
#include <fstream> 
#include <string>

#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"

void compare_eDIS_eAll(string filename1,string filename2){
   gStyle->SetPalette(1);
   gStyle->SetOptStat(0);

   TF1 *W_2=new TF1("W_2","(0.938*0.938-4+2*0.938*11)/(2*0.938+4*11*sin(x/180.*3.1415926/2)*sin(x/180.*3.1415926/2))",0,200);
   W_2->SetLineColor(kRed);
   TF1 *Q2_1=new TF1("Q2_1","1./(4*11*sin(x/180.*3.1415926/2)*sin(x/180.*3.1415926/2))",0,200);  
   Q2_1->SetLineColor(kBlack);

   //TFile *file_eicrate=new TFile("/work/halla/solid/zwzhao/solid/solid_svn/solid/evgen/eicRate_20101102/output/rate_solid_PVDIS_LD2_eDIS_1e6.root");
   // TFile *file_eicrate=new TFile("/work/halla/solid/zwzhao/solid/solid_svn/solid/evgen/eicRate_20101102/output/rate_solid_PVDIS_LD2_eDIS_50GeV_1e6.root");
   TFile *file_eicrate=new TFile(filename1.c_str());

   TTree *T_eicrate = (TTree*) file_eicrate->Get("T");

   // TFile *file_solid_inclusive_e=new TFile("/work/halla/solid/zwzhao/solid/solid_svn/solid/evgen/solid_inclusive_e/gen_solid_PVDIS_LD2_mode1_1e6.root");
   //TFile *file_solid_inclusive_e=new TFile("/work/halla/solid/zwzhao/solid/solid_svn/solid/evgen/solid_inclusive_e/gen_solid_PVDIS_LD2_mode1_1e6_cteq66.root");
   TFile *file_solid_inclusive_e=new TFile(filename2.c_str());

   // TFile *file_solid_inclusive_e=new TFile("/work/halla/solid/zwzhao/solid/solid_svn/solid/evgen/solid_inclusive_e/gen_solid_PVDIS_LD2_mode1_50GeV_1e6.root");

   TTree *T_solid_inclusive_e = (TTree*) file_solid_inclusive_e->Get("T");

   TCanvas *c_rate = new TCanvas("compare_this","compare_this",1800,900);
   c_rate->Divide(3,1);

   c_rate->cd(1);
   gPad->SetLogz();
   // T_eicrate->Draw("pf:theta/3.1415926*180>>heicrate(180,0,180,110,0,11)","rate","colz");
   // heicrate->SetTitle("eDIS rate;theta(deg);mom(GeV)");
   // W_2->Draw("same");
   // Q2_1->Draw("same");
   T_eicrate->Draw("Q2:W>>heicrate(100,0,5,100,0,20)","rate*(W<3&&x>0.35&&theta/3.1415926*180>22&&theta/3.1415926*180<35)","colz"); 
   // T_eicrate->Draw("Q2:W>>heicrate(120,0,12,100,0,100)","rate*(W>3&&Q2>1&&Q2>0.2&&theta/3.1415926*180<50)","colz"); 
   heicrate->SetTitle("eDIS rate;W(GeV);Q2(GeV^{2})");
   // T_eicrate->Draw("Q2:W>>heicrate(120,0,12,100,0,100)","rate*(W>3&&Q2>1&&Q2>0.2&&theta/3.1415926*180<50)","colz"); 
   // heicrate->SetTitle("eDIS rate;W(GeV);Q2(GeV^{2})");
   // T_eicrate->Draw("Q2:x>>heicrate(100,0,1,100,0,100)","rate*(W>3&&Q2>1&&Q2>0.2&&theta/3.1415926*180<50)","colz"); 
   // heicrate->SetTitle("eDIS rate;x;Q2(GeV^{2})");
   heicrate->SetMaximum(1e13);
   heicrate->SetMinimum(1e-3);

   c_rate->cd(2);
   gPad->SetLogz();
   // T_solid_inclusive_e->Draw("sqrt(px*px+py*py+pz*pz):theta>>h_solid_inclusive_e(180,0,180,110,0,11)","rate","colz");
   // h_solid_inclusive_e->SetTitle("eAll rate;theta(deg);mom(GeV)");
   // W_2->Draw("same");
   // Q2_1->Draw("same");
   T_solid_inclusive_e->Draw("Q2:W>>h_solid_inclusive_e(100,0,5,100,0,20)","rate*(W>3&&x>0.35&&theta>22&&theta<35)","colz");
   // T_solid_inclusive_e->Draw("Q2:W>>h_solid_inclusive_e(120,0,12,100,0,100)","rate*(W>3&&Q2>1&&Q2>0.2&&theta<50)","colz");
   h_solid_inclusive_e->SetTitle("eAll rate;W(GeV);Q2(GeV^{2})");
   // T_solid_inclusive_e->Draw("Q2:x>>h_solid_inclusive_e(100,0,1,100,0,100)","rate*(W>3&&Q2>1&&Q2>0.2&&theta<50)","colz");
   // h_solid_inclusive_e->SetTitle("eAll rate;x;Q2(GeV^{2})");
   h_solid_inclusive_e->Draw("colz");
   h_solid_inclusive_e->SetMaximum(1e13);
   h_solid_inclusive_e->SetMinimum(1e-3);


   c_rate->cd(3);
   gPad->SetLogz();
   TH2F *hcompare=heicrate->Clone();
   // hcompare->Divide(heicrate,h_solid_inclusive_e);
   // hcompare->SetTitle("eDIS/eAll");
   hcompare->Divide(h_solid_inclusive_e,heicrate);
   hcompare->SetTitle("eAll/eDIS");
   hcompare->SetMaximum(100);
   hcompare->SetMinimum(0);
   hcompare->Draw("colz");
   // W_2->Draw("same");
   // Q2_1->Draw("same");

   // 22-35 deg
   // cout << "eDIS rate "<< heicrate->Integral(44,70,10,110) << endl;
   // cout << "eAll rate "<< h_solid_inclusive_e->Integral(44,70,10,110) << endl;

   cout << "eDIS rate "<< heicrate->Integral() << endl;
   cout << "eAll rate "<< h_solid_inclusive_e->Integral() << endl;

}
