#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <THnSparse.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>
#include <TVector3.h>
#include "TRandom3.h"
#include "TSystem.h"
#include "TLorentzRotation.h"

using namespace std;

// void analysis_eAll(string input_filename,string acc_filename)
void analysis_eAll(string input_filename,string acc_filename,double beam)
{
gROOT->Reset();
// gStyle->SetPalette(1);
gStyle->SetOptStat("i");
gStyle->SetOptFit(0);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
//   gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadLeftMargin(0.2);  
  gStyle->SetPadRightMargin(0.15);  
  
  gStyle->SetPadColor(0);

  gStyle->SetLabelSize(0.04,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.04,"xyz");   
  gStyle->SetTitleSize(0.07,"t");    
//   gStyle->SetPaintTextFormat("4.1f"); 

const double DEG=180./3.1415926;

// char beam_s[200];
// int beam;
// if(sprintf(beam_s,"%s",input_filename.substr(input_filename.find("GeV")-2,2).c_str())) {
// // cout << beam_s << endl;
// beam=atoi(beam_s);
// cout << "beam energy " << beam << endl;
// }
// else {cout << "wrong beam energy" << endl; return;}

bool Is_NH3=false;
if (input_filename.find("NH3",0) != string::npos) {
    Is_NH3=true;
    cout << "SIDIS NH3 using 3D acceptance" << endl;
}

TFile *file_acc = new TFile(acc_filename.c_str());  
TH2F *hf,*hl,*hall;
hf = (TH2F*)file_acc->Get("acceptance_ThetaP_forwardangle");
hl = (TH2F*)file_acc->Get("acceptance_ThetaP_largeangle");
hall = (TH2F*)file_acc->Get("acceptance_ThetaP_overall");
TH3F *hf_3D,*hl_3D;
hf_3D = (TH3F*)file_acc->Get("acceptance_ThetaPhiP_forwardangle");
hl_3D = (TH3F*)file_acc->Get("acceptance_ThetaPhiP_largeangle");

TCanvas *c_acc = new TCanvas("c_acc","c_acc",1800,900);
c_acc->Divide(2,1);
c_acc->cd(1);
gPad->SetLogz();
if(!Is_NH3) hf->Draw("colz");
else hf_3D->Draw("box");
c_acc->cd(2);
gPad->SetLogz();
if(!Is_NH3) hl->Draw("colz");
else hl_3D->Draw("box");

TH2F *hrate_Theta_P_gen = new TH2F("hrate_Theta_P_gen","rate (kHz);vertex #theta (deg);vertex P (GeV)",50,0,50,int(beam*10),0,beam);
TH2F *hrate_x_Q2_gen = new TH2F("hrate_x_Q2_gen","rate (kHz);x_{bj};Q^{2} (GeV^{2});",100,0,1,int(beam*20),0,beam*2);

TH2F *hxs_Theta_P_gen = new TH2F("hxs_Theta_P_gen","crosssection (#mub);vertex #theta (deg);vertex P (GeV)",50,0,50,int(beam*10),0,beam);
TH2F *hxs_x_Q2_gen = new TH2F("hxs_x_Q2_gen","crosssection (#mub);x_{bj};Q^{2} (GeV^{2});",100,0,1,int(beam*20),0,beam*2);

TH2F *hrate_Theta_P = new TH2F("hrate_Theta_P","rate (kHz);vertex #theta (deg);vertex P (GeV)",50,0,50,int(beam*10),0,beam);
TH2F *hrate_x_Q2 = new TH2F("hrate_x_Q2","rate (kHz);x_{bj};Q^{2} (GeV^{2});",100,0,1,int(beam*20),0,beam*2);

TH2F *hxs_Theta_P = new TH2F("hxs_Theta_P","crosssection (#mub);vertex #theta (deg);vertex P (GeV)",50,0,50,int(beam*10),0,beam);
TH2F *hxs_x_Q2 = new TH2F("hxs_x_Q2","crosssection (#mub);x_{bj};Q^{2} (GeV^{2});",100,0,1,int(beam*20),0,beam*2);

TFile *file=new TFile(input_filename.c_str());
if (file->IsZombie()) {
    cout << "Error opening file" << input_filename << endl;
    exit(-1);
}
else cout << "open file " << input_filename << endl;    

TTree *tree = (TTree*) file->Get("T");

   Double_t        Abeam;
   Double_t        AL;
   Double_t        x;
   Double_t        y;
   Double_t        W;
   Double_t        Q2;
   Int_t           charge;
   Int_t           particle_id;
   Double_t        px;
   Double_t        py;
   Double_t        pz;
   Double_t        Ep;
   Double_t        mass;
   Double_t        vx;
   Double_t        vy;
   Double_t        vz;
   Double_t        xs;
   Double_t        theta;
   Double_t        phi;
   Double_t        rate;

   // List of branches
   TBranch        *b_data;
   tree->SetBranchAddress("Abeam", &Abeam, &b_data);
   tree->SetBranchAddress("AL", &AL, &b_data);
   tree->SetBranchAddress("x", &x, &b_data);
   tree->SetBranchAddress("y", &y, &b_data);
   tree->SetBranchAddress("W", &W, &b_data);
   tree->SetBranchAddress("Q2", &Q2, &b_data);
   tree->SetBranchAddress("charge", &charge, &b_data);
   tree->SetBranchAddress("particle_id", &particle_id, &b_data);
   tree->SetBranchAddress("px", &px, &b_data);
   tree->SetBranchAddress("py", &py, &b_data);
   tree->SetBranchAddress("pz", &pz, &b_data);
   tree->SetBranchAddress("Ep", &Ep, &b_data);
   tree->SetBranchAddress("mass", &mass, &b_data);
   tree->SetBranchAddress("vx", &vx, &b_data);
   tree->SetBranchAddress("vy", &vy, &b_data);
   tree->SetBranchAddress("vz", &vz, &b_data);
   tree->SetBranchAddress("xs", &xs, &b_data);
   tree->SetBranchAddress("theta", &theta, &b_data);
   tree->SetBranchAddress("phi", &phi, &b_data);
   tree->SetBranchAddress("rate", &rate, &b_data);
   
int nevent = (int)tree->GetEntries();
// double time = 48.0 * 24.0 * 3600.0;
cout << "nevent " << nevent << endl;
for (Int_t i=0;i<nevent;i++) { 
//   cout << i << "\r";
  tree->GetEntry(i);
       double p=sqrt(px*px+py*py+pz*pz);
//        cout << Ep << " " << p << endl;  //they are same
       double theta=acos(pz/Ep)*DEG;
       double phi=atan2(px,py)*DEG;       
//        if(phi<-180 || phi>180) cout<<phi<<endl; 

       int bin,binx,biny,binz;
       double acc_f, acc_l, acc;
       if(!Is_NH3){
        bin=hf->FindBin(theta,p);
        hf->GetBinXYZ(bin,binx,biny,binz);
        acc_f = hf->GetBinContent(binx,biny); 
        bin=hl->FindBin(theta,p);
        hl->GetBinXYZ(bin,binx,biny,binz);
        acc_l = hl->GetBinContent(binx,biny); 
       }
       else{
        bin=hf->FindBin(theta,phi,p);
        hf->GetBinXYZ(bin,binx,biny,binz);
        acc_f = hf->GetBinContent(binx,biny,binz); 
        bin=hl->FindBin(theta,phi,p);
        hl->GetBinXYZ(bin,binx,biny,binz);
        acc_l = hl->GetBinContent(binx,biny,binz); 
       }
       
//        apply some cut if needed
//       if(Q2>1) acc_f=0;       
//       if(p>3.0) acc_l=0;
//       if(p>3.0 && W>2 && Q2>9 && Q2<10){}
//       else acc_l=0;
//       rate=rate*1e3*3600*24*48/1e8;
       
//       acc=acc_l;
      acc=acc_f+acc_l;
      
      hrate_Theta_P_gen->Fill(theta,p,rate/1e3);       //in khz
      hrate_x_Q2_gen->Fill(x,Q2,rate/1e3);       //in khz
      hxs_Theta_P_gen->Fill(theta,p,xs);       //in ub
      hxs_x_Q2_gen->Fill(x,Q2,xs);       //in ub

      hrate_Theta_P->Fill(theta,p,rate*acc/1e3);       //in khz
      hrate_x_Q2->Fill(x,Q2,rate*acc/1e3);       //in khz
      hxs_Theta_P->Fill(theta,p,xs*acc);       //in ub
      hxs_x_Q2->Fill(x,Q2,xs*acc);       //in ub

}
file->Close();

TCanvas *c_rate_gen = new TCanvas("c_rate_gen","c_rate_gen",1800,900);
c_rate_gen->Divide(2,1);
c_rate_gen->cd(1);
gPad->SetLogz();
hrate_Theta_P_gen->Draw("colz");
c_rate_gen->cd(2);
// gPad->SetLogz();
hrate_x_Q2_gen->Draw("colz");
c_rate_gen->SaveAs("rate_gen.png");
cout << "total rate gen " << hrate_Theta_P_gen->Integral() << " khz" << endl;

TCanvas *c_xs_gen = new TCanvas("c_xs_gen","c_xs_gen",1800,900);
c_xs_gen->Divide(2,1);
c_xs_gen->cd(1);
gPad->SetLogz();
hxs_Theta_P_gen->Draw("colz");
c_xs_gen->cd(2);
// gPad->SetLogz();
hxs_x_Q2_gen->Draw("colz");
c_xs_gen->SaveAs("xs_gen.png");
cout << "total xs gen " << hxs_Theta_P_gen->Integral() << " ub" << endl;

TCanvas *c_rate = new TCanvas("c_rate","c_rate",1800,900);
c_rate->Divide(2,1);
c_rate->cd(1);
gPad->SetLogz();
// hrate_Theta_P->SetMinimum(0.1);
hrate_Theta_P->Draw("colz");
c_rate->cd(2);
// gPad->SetLogz();
// hrate_x_Q2->SetMinimum(0.1);
hrate_x_Q2->Draw("colz");
c_rate->SaveAs("rate_acc.png");
cout << "total rate " << hrate_Theta_P->Integral() << " khz" << endl;

TCanvas *c_xs = new TCanvas("c_xs","c_xs",1800,900);
c_xs->Divide(2,1);
c_xs->cd(1);
gPad->SetLogz();
// hxs_Theta_P->SetMinimum(0.1);
hxs_Theta_P->Draw("colz");
c_xs->cd(2);
// gPad->SetLogz();
// hxs_x_Q2->SetMinimum(0.1);
hxs_x_Q2->Draw("colz");
c_xs->SaveAs("xs_acc.png");
cout << "total xs " << hxs_Theta_P->Integral() << " ub" << endl;

TFile *outputfile_final=new TFile("eAll_output.root", "recreate");
outputfile_final->Append(hrate_Theta_P_gen);
outputfile_final->Append(hrate_x_Q2_gen);
outputfile_final->Append(hxs_Theta_P_gen);
outputfile_final->Append(hxs_x_Q2_gen);
outputfile_final->Append(hrate_Theta_P);
outputfile_final->Append(hrate_x_Q2);
outputfile_final->Append(hxs_Theta_P);
outputfile_final->Append(hxs_x_Q2);
outputfile_final->Write();
outputfile_final->Flush();

}
