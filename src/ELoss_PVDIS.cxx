//g2p energy loss
#include <algorithm>
#include <typeinfo>
#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <string.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TChain.h>
#include <TMath.h>
#include <TBenchmark.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TPaveStats.h>
#include <TImage.h>

#include "ELoss_PVDIS.h"


////////////////////////////////////////////////////////
namespace PVDIS {

    //pTargetType: 1 is LH2 2 is LD2
    void ConstructTarget(int pTgtType, double pTgtEntrZPos)
    {
        //PVDIS target contains the folowing"
        //version 1: beam pipe does not connect to the target chamber
        //1) 10-mil Be beam exit, 
        //2) then 10cm air or helium4 gas, 
        //3) then 10-mil aluminum target chamber entrance window, 
        //4) then 4-mil aluminum upstream target cap, 
        //5) then 40 cm LH2 or LD2 target, 
        //6) then 5 mil aluminum downstream target cap

        //version 2: beam pipe connect to the target chamber
        //1) 10-mil Be beam exit, 
        //2) then 4-mil aluminum upstream target cap, 
        //3) then 40 cm LH2 or LD2 target, 
        //4) then 5 mil aluminum downstream target cap

        //Here I used version 1, but item 6) will not build

        gTgtType=pTgtType;
        gTgtEntrZPos=pTgtEntrZPos;

        const double ZH=1; //Z
        const double ZD=1; //Z
        const double ZHe=2;
        const double ZBe=4;
        const double ZN=7;
        const double ZO=8;
        const double ZAl=13;
        const double ZAir=7;

        const double AH=1.0087; //A
        const double AD=2.014101764; //A
        const double AHe=4.002602;
        const double ABe=9.012182;
        const double AN=14.00674;
        const double AO=15.9994;
        const double AAl=26.981539;
        const double AAir=14.421;

        //density,g/cm3
        double DLH2=0.0708 * g/cm3;
        double DLD2=0.1638 * g/cm3;
        double DHe=0.000164 * g/cm3;   //gas in 1 ATM     
        double DLHe=0.145 * g/cm3;
        double DBe=1.848 * g/cm3;
        double DLN2=0.807 * g/cm3;
        double DBeO=3.02 * g/cm3;
        double DAir=0.000129 * g/cm3;
        double DAl=2.6989 * g/cm3; 

        //Rad.len. g/cm2 X0, 
        double XH=63.04 * g/cm2;
        double XD=125.98 * g/cm2;
        double XHe=94.322 * g/cm2; 
        double XBe=65.19 * g/cm2;
        double XN=40.7 * g/cm2;
        double XAl=24.011 * g/cm2;
        double XAir=36.66 * g/cm2;
        
	//thickness
        double LBeamExit=10 * mil;
        double LAirEntr=10. * cm;
        double LTgtChamEntr=10 * mil;
        double LTgtEntrCap=4 *mil;
        double LTgtLengh=40 * cm;
        double LTgtExitCap=5 *mil;


        int idx=0;

        //BeO beam exit
        //int nCompnent=2;
        //double Zlist_BeO[]={ZBe,ZO};
        //double Alist_BeO[]={ABe,AO};
        //int Nlist_BeO[]={1,1};
        //ConstructMaterial(gMatl[idx++],nCompnent,Zlist_BeO,Alist_BeO,Nlist_BeO,DBeO,LBeamExit,"BeO","BeamExit");
        //or Be beam exit??
        ConstructMaterial(gMatl[idx++],ZBe,ABe,DBe,LBeamExit,"Be","BeamExit",0,XBe);

        //Helium gas in 1 ATM
        //ConstructMaterial(gMatl[idx++],ZHe,AHe,DHe,LTgtNH3,"He4Gas","He4Gas",0,XHe);
        //Air entrance, 1 ATM
        ConstructMaterial(gMatl[idx++],ZAir,AAir,DAir,LAirEntr,"Air","AirEntr",0,XAir);


        //target chamber entrance
        ConstructMaterial(gMatl[idx++],ZAl,AAl,DAl,LTgtChamEntr,"Al","TgtChamEntr",0,XAl);

        //Upstream endcap 
        ConstructMaterial(gMatl[idx++],ZAl,AAl,DAl,LTgtEntrCap,"Al","TgtEntrCap",0,XAl);

        //LH2 or LD2
        if(gTgtType==1) {
            ConstructMaterial(gMatl[idx++],ZH,AH,DLH2,LTgtLengh,"LH2","TgtLH2",0,XH); 
        } 
        else {
            ConstructMaterial(gMatl[idx++],ZD,AD,DLD2,LTgtLengh,"LD2","TgtLD2",0,XD);
        }

        //target exit cap (downstream end cap)
        //ConstructMaterial(gMatl[idx++],ZAl,AAl,DAl,LTgtExitCap,"Al","TgtExitCap",0,XAl);

        //LHe4, not use here
        //ConstructMaterial(gMatl[idx++],ZHe,AHe,DLHe,LLHe,"He","TgtHe",IHe,XHe);

        //Vacuum
        //ConstructMaterial(gMatl[idx++],0,0,0,0,"Vacuum","Vacuum",0,999999999.0);

        SetRandomSeed(0);

        cout<<"\n PVDIS "<<((gTgtType==1)?"LH2":"LD2")<<" target is constructed:\n";
        cout<<setw(15)<<"GeoName"<<" "<<setw(8)<<"Material"<<" "<<setw(8)<<"Thick(cm)"<<" "<<setw(13)<<"DL/X0"<<" "<<setw(8)<<"D(g/cm3)"<<endl;
        for(int i=0;i<kNMatl;i++)
        {
            cout<<setw(15)<<gMatl[i].name<<" "<<setw(8)<<gMatl[i].matl_name
                <<" "<<setw(8)<<gMatl[i].L;
            cout.setf(std::ios::scientific);    
            cout<<" "<<setw(13)<<gMatl[i].DL/gMatl[i].X0;
            cout.unsetf(std::ios::scientific);
            cout<<" "<<setw(8)<<gMatl[i].D
                <<endl;
        }
    }

    //calculate energy after Mean energy loss for proton
    double CalELoss_proton(double E0, double vz, double *Ef)
    {
        //change the thickness of target according to given vz
        gMatl[kNMatl-1].L = vz - gTgtEntrZPos;
        
        for(int i=0;i<kNMatl;i++) 
        {
            Ef[i]=(i>0)?Ef[i-1]:E0;
            if(gMatl[i].A>0.0 && Ef[i]>0.0)
            {
                Ef[i]=ELOSS::GetNewE_heavy(Ef[i],gMatl[i],M_p,1);
            }
        }

        return Ef[kNMatl-1];
    }

    //calculate energy after energy loss for electron
    double CalELoss_electron(double E0, double vz, int WithIntRad, double *Ef)
    {
        //change the thickness of target according to given vz
        gMatl[kNMatl-1].L = vz - gTgtEntrZPos;
        double Angle=0.0;  //for beam electron, ignoring MSC, therefore it is zero
        
        for(int i=0;i<kNMatl;i++)
        {
            Ef[i]=(i>0)?Ef[i-1]:E0;
            if(gMatl[i].A>0.0 && Ef[i]>0.0)
            {
                Ef[i]=ELOSS::GetNewE(Ef[i],gMatl[i],WithIntRad, Angle);
            }
        }

        return Ef[kNMatl-1];
    }

    //calculate energy after energy loss for electron
    double CalELoss_electron(double E0, double vz, int WithIntRad)
    {       
       double Ef[kNMatl];
       return CalELoss_electron( E0, vz, WithIntRad, Ef);
    }

    //plot energy loss figure
    void PlotELoss(double beamE, int nThrown=1000)
    {
        int i,j;
        static const int kNBins=200; //number of energy bins
        const double kEBinWidth=(beamE+0.1)/double(kNBins);
        double **Eall=0; //Eall[nThrown][kNMatl]
        Eall = new double * [nThrown];
        for(i=0;i<nThrown;i++)
        {
            Eall[i] = new double [kNMatl];
            for(j=0;j<kNMatl;j++) Eall[i][j]=0.0;
        }

        double E_life_ave[kNMatl],E_live_sum[kNMatl];
        double N_E[kNMatl][kNBins];//number of electrons in every energy bin
        double E_E[kNMatl][kNBins]; //energy in each energy bin
        double N_live[kNMatl];//number of electrons pass through material

        system("mkdir Graph");
        ofstream EallFile("Graph/Eall.txt");
	
        EallFile<<"Event     Ef[0]     Ef[1]     Ef[2]     Ef[3]     Ef[4]\n";
        int pWithIntRad=0;
        for(i=0;i<kNMatl;i++)
        {
            E_live_sum[i]=N_live[i]=E_life_ave[i]=0.0;
            //initialize the arrays which will be used to plot graphs
            for(j=0;j<kNBins;j++)
            {
                N_E[i][j]=0;
                E_E[i][j]=kEBinWidth*(j+0.5);
            }
        }

        double Ee[kNMatl];
        for(i=0;i<nThrown;i++)
        {
            Ee[0]=beamE;
            for (j=1;j<kNMatl;j++) Ee[j]=0;  //reset
            double vz=10*cm;

            //calculate energy after energy loss for electron
            //double CalELoss_electron(double E0, double vz, int WithIntRad, double *Ef)
            CalELoss_electron(beamE, vz, pWithIntRad, Ee);

            EallFile<<setw(5)<<i<<"  ";
            for(j=0;j<kNMatl;j++)
            {
                Eall[i][j]=Ee[j];
                if(Ee[j]>0.0)
                {
                    E_live_sum[j]+=Ee[j];
                    N_live[j]++;
                    //calculate the number of electrons in each energy bin
                    int bin=int(Ee[j]/kEBinWidth);
                    if(bin>=0 && bin<kNBins) N_E[j][bin]++;
                }
                else continue;
                EallFile<<setw(8)<<Eall[i][j]<<"  ";
            }
            EallFile<<endl;
        }
        EallFile.close();

        double Rx[kNMatl];
        Rx[0]=gMatl[0].L;
        for(int i=1;i<kNMatl;i++)
        {
            Rx[i]=Rx[i-1]+gMatl[i].L;
        }

        ofstream EaveFile("Graph/E_life_ave.txt");
        EaveFile<<"       Material   D*L/X0    E(MeV) TrackL(cm)  Count\n";
        cout    <<"       Material   D*L/X0    E(MeV) TrackL(cm)  Count\n";

        char tmpStr[200];
        sprintf(tmpStr,"%15s  %7.4f  %8.3f  %8.3f  %6d\n","Entrance",0.0,beamE,0.0,nThrown);
        EaveFile<<tmpStr;
        cout<<tmpStr;

        for(i=0;i<kNMatl;i++)
        {
            E_life_ave[i]=(N_live[i]>0.0)?E_live_sum[i]/N_live[i]:0.0;
            sprintf(tmpStr,"%15s  %7.4f  %8.3f  %8.3f  %6d\n",gMatl[i].name.c_str(),
                gMatl[i].DL/gMatl[i].X0,E_life_ave[i],Rx[i],int(N_live[i]));
            EaveFile<<tmpStr;
            cout<<tmpStr;
        }

        EaveFile.close();

        /////////////////////// mean energy loss ////////////////////////////
        TCanvas *c1=new TCanvas("c1","Mean Energy Loss",800,600);
        TGraph *g1=new TGraph(kNMatl,Rx,E_life_ave);
        g1->GetXaxis()->SetTitle("Position(cm)");
        g1->GetYaxis()->SetTitle("Energy(MeV)");
        g1->SetTitle("Mean Energy Loss");
        g1->Draw("AL*");

        c1->SaveAs("Graph/E_life_ave.png");
        c1->SaveAs("Graph/E_life_ave.eps");


        /////////////////////// energy loss distribution ///////////////////////
        TGraph *gr[kNMatl];
        TAxis *axisX[kNMatl];
        TCanvas *c2=new TCanvas("c2","distribution",1280,900);
        c2->Divide(3,2,0.001,0.001);

        for(j=0;j<kNMatl;j++)
        {
            c2->cd(j+1);
            gr[j]=new TGraph(kNBins,E_E[j],N_E[j]);
            gr[j]->GetXaxis()->SetTitle("Energy(MeV)");
            gr[j]->GetYaxis()->SetTitle("Number of Electrons");
            sprintf(tmpStr,"%s(DL/X=%5.3f,I=%5.1feV):%4d in,%4d out",gMatl[j].name.c_str(),
                gMatl[j].DL/gMatl[j].X0,gMatl[j].I/eV,int((j>0)?N_live[j-1]:nThrown),int(N_live[j]));
            gr[j]->SetTitle(tmpStr);
            gr[j]->Draw("AL");
	    gPad->SetLogy();
	}
        axisX[0]=gr[0]->GetXaxis();
        axisX[0]->SetLimits(beamE*2./3.,beamE*4./3.);
        axisX[1]=gr[1]->GetXaxis();
        axisX[1]->SetLimits(beamE*2./3.,beamE*4./3.);
        axisX[2]=gr[2]->GetXaxis();
        axisX[2]->SetLimits(beamE*2./3.,beamE);
        axisX[3]=gr[3]->GetXaxis();
        axisX[3]->SetLimits(beamE*2./3.,beamE);
        axisX[4]=gr[4]->GetXaxis();
        axisX[4]->SetLimits(beamE*2./3.,beamE);

        gr[0]->SetLineColor(1);
        gr[1]->SetLineColor(3);
        gr[2]->SetLineColor(4);
        gr[3]->SetLineColor(5);
        gr[4]->SetLineColor(2);

        gr[0]->SetMarkerColor(1);
        gr[1]->SetMarkerColor(3);
        gr[2]->SetMarkerColor(4);
        gr[3]->SetMarkerColor(5);
        gr[4]->SetMarkerColor(2);

        gr[0]->SetLineWidth(3);
        gr[1]->SetLineWidth(3);
        gr[2]->SetLineWidth(3);
        gr[3]->SetLineWidth(3);
        gr[4]->SetLineWidth(3);

        c2->SaveAs("Graph/Eall.png");
        //c2->SaveAs("Graph/Eall.eps");
        //c2->SaveAs("Graph/Eall.c");

        //////////////////////////////////////////////////////////////
        TCanvas *c3=new TCanvas("c3","Energy Loss Distribution",800,600);
        TMultiGraph *md=new TMultiGraph();
        TLegend* leg1 = new TLegend(0.14,0.5,0.45,0.89);
	leg1->SetFillColor(0);
	if(gROOT->IsBatch()) leg1->SetFillStyle(4000);
        const int kNPlot=5;
        int PlotGrIdx[]={0,1,2,3,4};

        for(j=0;j<kNPlot;j++)
        {
            int idx=PlotGrIdx[j];
            gr[idx]->SetLineWidth(3);
            md->Add(gr[idx],"c");
            sprintf(tmpStr,"%s: %4d",gMatl[idx].name.c_str(),int(N_live[idx]));
            leg1->AddEntry(gr[idx],tmpStr,"l");
        }

        md->SetTitle("Distribution of Energy");
	c3->cd();
        c3->SetLogy();
        md->Draw("A");
        leg1->Draw("same");
        md->GetXaxis()->SetTitle("Energy(MeV)");
        md->GetYaxis()->SetTitle("Number of electrons");
	
        c3->SaveAs("Graph/Edistr.png");
        //c3->SaveAs("Graph/Edistr.C");


        //free the memory
        for(i=0;i<nThrown;i++) delete [] Eall[i];

        delete g1;
        for(j=0;j<kNMatl;j++) delete gr[j];
        delete md;
        delete c1;
        delete c2;
        delete c3;
    }

};  //end of namespace PVDIS