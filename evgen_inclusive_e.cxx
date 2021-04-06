#include<iostream>
#include<fstream>
#include<cmath>
#include "math.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"

//defined in the include directory
#include "constants.h"
#include "christy_bosted_inelastic_QE.h"
#include "proton_DIS.h"
#include "neutron_DIS.h"
#include "fixed_target_xs.h"
#include "norad_XS.h"
#include "RadiativeCorrections.h"
#include "eInclusiveCrossSection.h"
#include "InputParameters.h"
#include "FileManager.h"

#include "ELoss_PVDIS.h"

using namespace LHAPDF;
using namespace std;

//############################################################################################
//             Summary of functions can be used
/*
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

//proton structure functions
double calculate_proton_g1gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g5gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g3gz(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g5z(PDF* pol_pdf, double x, double Q2);

double calculate_proton_g3z(PDF* pol_pdf, double x, double Q2);

double calculate_proton_F2g(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F1g(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F1gz(PDF* unpol_pdf, double x, double Q2);

double calculate_proton_F3gz(PDF* unpol_pdf, double x, double Q2);

//---proton asymmetries
//---PVDIS asymmetry using unpol electron + longitudinally polarized proton
double calculate_proton_AL(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);

double calculate_proton_AL_g1gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);

double calculate_proton_AL_g5gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);

//---PVDIS asymmetry using long. pol. electron + unpolarized proton
double calculate_proton_Abeam(PDF* unpol_pdf, double x, double Q2, double y);


double calculate_fixed_target_xs(double E, int Z, int A, double theta, double Ep, PDF* unpol_pdf);

//Peter Bostd model
int F1F2IN09(int Z, int IA, double qsq, double wsq, double &F1, double &F2, double &Rc);

// Peter Bosted model
void F1F2QE09(int Z, int IA, double QSQ, double wsq, double &F1, double &F2);

*/
//############################################################################################


int  main(Int_t argc, char *argv[])
{
    if(argc<2) {
        cout<<"***Error: You need to provide at least 1 argument\n";
        cout<<"Usage: "<<argv[0]<<" <input_file> [doIncomingEloss=1] [doRadiateBornXS=1]\n";
        exit(-1);
    }
    char input_gen_file[250];
    strcpy(input_gen_file,argv[1]);
    std::string str = input_gen_file;
    std::cout<<"input file "<< str <<std::endl;
    std::cout<<std::endl;

    int doIncomingEloss = 1;
    if(argc>2) doIncomingEloss = atoi(argv[2]);
    int doRadiateBornXS = 1;
    if(argc>3) doRadiateBornXS = atoi(argv[3]);


    inputParameters_t par;
    FileManager *myFN = new FileManager();
    myFN->LoadInputData(input_gen_file,par);

    std::cout << "INPUT PARAMETERS: " << std::endl;
    std::cout << "num_evt:        " <<  par.num_evt              << std::endl;
    std::cout << "output_name:    " <<  par.output_name          << std::endl;
    std::cout << "pol PDF name:   " <<  par.pol_pdfset_name      << std::endl;
    std::cout << "pol PDF ID:     " <<  par.pol_pdfset_ID        << std::endl;
    std::cout << "unpol PDF name: " <<  par.unpol_pdfset_name    << std::endl;
    std::cout << "unpol PDF ID:   " <<  par.unpol_pdfset_ID      << std::endl;
    std::cout << "tgt Z:          " <<  par.tgt_Z                << std::endl;
    std::cout << "tgt A:          " <<  par.tgt_A                << std::endl;
    std::cout << "luminosity:     " <<  par.lumi                 << std::endl;
    std::cout << "E_beam:         " <<  par.E_beam               << std::endl;
    std::cout << "th_min:         " <<  par.theta_min            << std::endl;
    std::cout << "th_max:         " <<  par.theta_max            << std::endl;
    std::cout << "Ep_min:         " <<  par.Ep_min               << std::endl;
    std::cout << "Ep_max:         " <<  par.Ep_max               << std::endl;
    std::cout << "vx_min:         " <<  par.vx_min               << std::endl;
    std::cout << "vx_max:         " <<  par.vx_max               << std::endl;
    std::cout << "vy_min:         " <<  par.vy_min               << std::endl;
    std::cout << "vy_max:         " <<  par.vy_max               << std::endl;
    std::cout << "vz_min:         " <<  par.vz_min               << std::endl;
    std::cout << "vz_max:         " <<  par.vz_max               << std::endl;
    std::cout << "scale:          " <<  par.scale                << std::endl;
    std::cout << "rad:            " <<  par.rad                  << std::endl;
    std::cout << "RL:            " <<  par.RL                  << std::endl;
    std::cout << "RLb:            " <<  par.RLb                  << std::endl;
    std::cout << "RLa:            " <<  par.RLa                  << std::endl;
    std::cout << "Fit_model:      " <<  par.Fit_model            << std::endl;

    delete myFN;

    //###################################################################################
    //
    //         user inputs for the simulation
    //
    //###################################################################################
    const double lumi=par.lumi;  //Hz/cm2  this must be luminosity per nucleus (corresponding to the "A" input), not per nucleon!
    const double E=par.E_beam;   //GeV, electron beam energy
    const int Z=par.tgt_Z;  //target nucleus Z
    const int A=par.tgt_A;  //target nucleus A   cross section is calculated per A
    const double theta_min=par.theta_min ;  //degree, define solid angle
    const double theta_max=par.theta_max;  //degree, define solid angle
    const double Ep_min=par.Ep_min;     //GeV, define momentum space
    const double Ep_max=par.Ep_max;    //GeV, define momentum space
    const double vertex_x_min=par.vx_min;  //cm, depend on the raster size
    const double vertex_x_max=par.vx_max;   //cm
    const double vertex_y_min=par.vy_min;  //cm
    const double vertex_y_max=par.vy_max;   //cm
    const double vertex_z_min=par.vz_min;   //cm
    const double vertex_z_max=par.vz_max;  //cm
    const int  num_evt=par.num_evt;    //number of event to generate within the phase space
    const int scale_status=par.scale;    //0:no scale factor on cross sections; 1: applying a scale factor on cross sections
    const int rad_status=par.rad;    //0 Born cross section only; 1:Radiative cross section with smearing; 2:Radiative cross section without smearing
    const double RL=par.RL;   //radiation length of target
    double RL_before=par.RLb;   //radiation length before target
    double RL_after=par.RLa;   //radiation length after target
    const int Fit_model=par.Fit_model;    //9---2009 Christy-Bosted Fit; 21-----2021 Christy's Fit 
    string pol_pdfset_name=par.pol_pdfset_name;    //pol. pdfset name
    string unpol_pdfset_name=par.unpol_pdfset_name;   // unpol. pdfset name
    //string unpol_pdfset_name="CT14lo";   // unpol. pdfset name
    //string unpol_pdfset_name="cteq66";   // unpol. pdfset name

    /*const double lumi=0.63e39;  //Hz/cm2  this must be luminosity per nucleus (corresponding to the "A" input), not per nucleon!
    const double E=4.74;   //GeV, electron beam energy
    const int Z=2;  //target nucleus Z
    const int A=3;  //target nucleus A   cross section is calculated per A
    const double theta_min=0 ;  //degree, define solid angle
    const double theta_max=180;  //degree, define solid angle
    const double Ep_min=0;     //GeV, define momentum space
    const double Ep_max=4.74;    //GeV, define momentum space
    const double vertex_x_min=-0.25;  //cm, depend on the raster size
    const double vertex_x_max=0.25;   //cm
    const double vertex_y_min=-0.25;  //cm
    const double vertex_y_max=0.25;   //cm
    const double vertex_z_min=-10;   //cm
    const double vertex_z_max=30;  //cm
    const int num_evt=int(atof(argv[2]));    //number of event to generate within the phase space
    string pol_pdfset_name="NNPDFpol11_100";    //pol. pdfset name
    string unpol_pdfset_name="CT14nlo";*/   // unpol. pdfset name
    //string unpol_pdfset_name="CT14lo";   // unpol. pdfset name
    //string unpol_pdfset_name="cteq66";   // unpol. pdfset name
    TString name_rootfile_output=par.output_name;   // name of the rootfile to save output data

    //output lund file
    ofstream OUTPUT_lund;
    TString name_rootfile_output_s=name_rootfile_output;
    OUTPUT_lund.open(name_rootfile_output_s.ReplaceAll("root","lund").Data());
    if(!OUTPUT_lund){
        cout<<"error! can't open lund output!"<<endl;
    }

    //###################################################################################
    //
    //             W<3 GeV, use chrity and Bosted fit
    //             W>3 GeV, USE LHAPDF6 TO ACCESS PDFs
    //
    //   load unpolarized pdf set and polarized pdf set in LHAPDF6
    //   PDF set is selected by the name of the pdfs
    //
    //
    //
    //###################################################################################

    //pol pdf sets
    //string pol_pdfset_name="NNPDFpol11_100";    //pdfset name
    string pol_pdfset_id=par.pol_pdfset_ID;         //central file
    const int pol_pdf_file_id=boost::lexical_cast<int>(pol_pdfset_id);  //nothing, just translate "0000" to 0, that's it

    //unpol pdf sets
    //string unpol_pdfset_name="CT14nlo";   //pdfset name
    //string unpol_pdfset_name="CT14nnlo";   //pdfset name
    //string unpol_pdfset_name="MMHT2014nlo68cl";   //pdfset name
    //string unpol_pdfset_name="MMHT2014nnlo68cl";   //pdfset name
    //string unpol_pdfset_name="MMHT2014lo68cl";   //pdfset name
    string unpol_pdfset_id=par.unpol_pdfset_ID;     //central file
    const int unpol_pdf_file_id=boost::lexical_cast<int>(unpol_pdfset_id); //nothing, just translate "0000" to 0, that's it

    cout<<"##########################################"<<endl;
    cout<<"# YX : pol_pdfset_name= "<<pol_pdfset_name<<"   pol_pdf_file_id= "<<pol_pdf_file_id<<endl;
    cout<<"# YX : unpol_pdfset_name= "<<unpol_pdfset_name<<"   unpol_pdf_file_id= "<<unpol_pdf_file_id<<endl;

    //------------------------------
    //load pol pdf and pdf set
    //------------------------------
    LHAPDF::PDF* pol_pdf=LHAPDF::mkPDF(pol_pdfset_name, pol_pdf_file_id);
    cout<<"# YX : polarized PDF set is loaded......using "<<pol_pdfset_name<<endl;

    //PDF set to do uncertainty calculation
    LHAPDF::PDFSet *pol_pdf_set=new LHAPDF::PDFSet(pol_pdfset_name);
    vector<LHAPDF::PDF*> pol_pdfs=pol_pdf_set->mkPDFs(); //get pdfs, for uncetainty estimation
    // int n_member=pol_pdf_set->size();
    //cout<<"ERROR CL LEVEL: "<<pol_pdf_set->errorConfLevel()<<endl;


    // load unpol pdf
    LHAPDF::PDF* unpol_pdf=LHAPDF::mkPDF(unpol_pdfset_name, unpol_pdf_file_id);
    cout<<"# YX : unpolarized PDF set is loaded......using "<<unpol_pdfset_name<<endl;
    // scale factor
    cout<<"Scale factor is loaded......using "<<scale_status<<endl;
    // radiative rate
    cout<<"Radiative rate is loaded......using "<<rad_status<<endl;
    cout<<"Fit model is loaded......using Christy-Bosted Fit "<<Fit_model<<endl;
    // radiative effects before vertex
    // 	cout<<"Tb is loaded......using "<<RL_before<<endl;
    // 	cout<<"Ta is loaded......using "<<RL_after<<endl;


    //print some information
    vector<int> pol_pids=pol_pdf->flavors();
    for(int i=0; i<int(pol_pids.size());i++){
        cout<<"# YX  : (pol pids) "<<i<<"	"<<pol_pids[i]<<endl;
    }

    vector<int> unpol_pids=unpol_pdf->flavors();
    for(int i=0; i<int(unpol_pids.size());i++){
        cout<<"# YX  : (unpol pids) "<<i<<"	"<<unpol_pids[i]<<endl;
    }
    
    //##################################################################################
    //Load energy loss namespace
    PVDIS::ConstructTarget(A,vertex_z_min);
    //##################################################################################
    //
    //
    //
    //               Define outputs
    //
    //##################################################################################

    //TTree to save
    double Abeam=0, AL=0, x=0, y=0, W=0, Q2=0, rate=0,raterad=0;
    double Ei=0, xi=0, yi=0, Wi=0, Q2i=0;  //Jixie: after Eloss variables
    //rate_pre=0;
    int charge=-1, particle_id=11;
    double px=0, py=0, pz=0;
    double Ep=0;
    double mass=0.511/1000.0;  //GeV
    double vx=0, vy=0, vz=0;
    double xs=0;
    double dXSdEdOmega_mubGeVSr=0; // differential cross section
    double radxs=0;
    double raddXSdEdOmega_mubGeVSr=0; // differential cross section with radiative correction
    double theta=0;
    double phi=0;
    double nu=0;
    double factor=0;
    TFile *myfile=new TFile(name_rootfile_output,"RECREATE");
    TTree *T=new TTree("T","T");
    T->Branch("x", &x, "x/D");
    T->Branch("y", &y, "y/D");
    T->Branch("W", &W, "W/D");
    T->Branch("Q2", &Q2, "Q2/D");
    //T->Branch("rate_pre", &rate_pre, "rate_pre/D");       //before normalized
    T->Branch("rate", &rate, "rate/D");
    T->Branch("raterad", &raterad, "raterad/D");
    T->Branch("Ei",&Ei, "Ei/D");     //energy right before scattering, after Eloss
    T->Branch("Abeam", &Abeam, "Abeam/D");
    T->Branch("AL", &AL, "AL/D");
    T->Branch("xi", &xi, "xi/D");    //Jixie: calculated using Ei instead of beam energy
    T->Branch("yi", &yi, "yi/D");
    T->Branch("Wi", &Wi, "Wi/D");
    T->Branch("Q2i", &Q2i, "Q2i/D");
    T->Branch("charge",&charge,"charge/I");
    T->Branch("particle_id",&particle_id,"particle_id/I");
    T->Branch("px",&px, "px/D");
    T->Branch("py",&py, "py/D");
    T->Branch("pz",&pz, "pz/D");
    T->Branch("Ep",&Ep, "Ep/D");
    T->Branch("mass",&mass, "mass/D");
    T->Branch("vx",&vx, "vx/D");
    T->Branch("vy",&vy, "vy/D");
    T->Branch("vz",&vz, "vz/D");
    T->Branch("xs",&xs, "xs/D");
    T->Branch("radxs",&radxs, "radxs/D");
    T->Branch("dXSdEdOmega_mubGeVSr",&dXSdEdOmega_mubGeVSr,"dXSdEdOmega_mubGeVSr/D");
    T->Branch("raddXSdEdOmega_mubGeVSr",&raddXSdEdOmega_mubGeVSr,"raddXSdEdOmega_mubGeVSr/D");
    T->Branch("theta",&theta, "theta/D");
    T->Branch("phi",&phi, "phi/D");
    T->Branch("nu" ,&nu, "nu/D");

    //###################################################################################
    //
    //
    //
    //
    //       start to throw events and calculate weight for each event
    //
    //
    //
    //###################################################################################
    TRandom3 rand;
    rand.SetSeed(0);

    // TBranch *brate = T->Branch("rate", &rate, "data/D");       //in unit Hz
    int count = 0;
    for(int i=0; i<num_evt; i++){
        //uniform vx, vy, vz
        vx=rand.Uniform(vertex_x_min, vertex_x_max);
        vy=rand.Uniform(vertex_y_min, vertex_y_max);
        vz=rand.Uniform(vertex_z_min, vertex_z_max);

        //---------phase space
        //uniform Ep
        Ep=rand.Uniform(Ep_min, Ep_max);
        //uniform phi, theta
        phi=rand.Uniform(0,360);  //0 degree to 360 degree
        theta=rand.Uniform(theta_min, theta_max);


        double d_omiga=2*PI*(cos(theta_min*deg_to_rad) - cos(theta_max*deg_to_rad));
        double d_E=Ep_max-Ep_min;

        px=Ep*sin(theta*deg_to_rad)*cos(phi*deg_to_rad);
        py=Ep*sin(theta*deg_to_rad)*sin(phi*deg_to_rad);
        pz=Ep*cos(theta*deg_to_rad);

        
        //caculate incoming energy loss, ture off internal brem.
        //repeat till the energy before scttering larger than the scattered energy 
        if(doIncomingEloss) {
            Ei = 0.0;
            int tryN = 0;
            while (Ei < Ep && (++tryN)<1000) 
            {
                Ei = PVDIS::CalELoss_electron(E*GeV, vz, 0)/1000.;  //turn into GeV
            }
        }
        else {
            Ei = E; 
        }

        //calculate after Eloss kinematics
        nu=Ei-Ep;
        Q2i=4.0*Ei*Ep*sin(theta*deg_to_rad/2.0)*sin(theta*deg_to_rad/2.0);
        Wi=sqrt(proton_mass*proton_mass + 2*proton_mass*nu - Q2i);
        xi=Q2i/2/proton_mass/nu;
        yi=nu/Ei;

        //calculate before Eloss kinematics
        nu=E-Ep;
        Q2=4.0*E*Ep*sin(theta*deg_to_rad/2.0)*sin(theta*deg_to_rad/2.0);
        W=sqrt(proton_mass*proton_mass + 2*proton_mass*nu - Q2);
        x=Q2/2/proton_mass/nu;
        y=nu/E;
        

        if(x>=0 && x<=1){
            //By Jixie: use the Ei to calculate XS now 
            //But still use beam to do the scaling since the scale factor is determined using the norminal beam
            xs=calculate_fixed_target_xs( Ei,  Z,  A,  theta,  Ep,  unpol_pdf, Fit_model);   //theta unit in degree
            if(scale_status==1){
                factor=0.906-0.00699*E;
            }else{
                factor=1.0;
            }
            //xs in unit of mub/GeV-sr
            //dXSdEdOmega_mubGeVSr = xs;
            //xs=xs*(d_E*d_omiga/num_evt);  //in unit of mub now
            eInclusiveCrossSection *noXS;
            noradXS NoRdXSDEdOmega;

            NoRdXSDEdOmega.SetZ(Z);
            NoRdXSDEdOmega.SetA(A);
            NoRdXSDEdOmega.SetEs(Ei);
            NoRdXSDEdOmega.SetEp(Ep);
            NoRdXSDEdOmega.SetTh(theta);
            NoRdXSDEdOmega.Setpdf(unpol_pdf);
            NoRdXSDEdOmega.SetScale(factor);
            // cout<<"dXSdEdOmega="<<dXSdEdOmega_mubGeVSr<<endl;
            noXS=&NoRdXSDEdOmega;
            double noradCross=noXS->GetBornXS();
            dXSdEdOmega_mubGeVSr = noradCross;
            // cout<<"noradCross="<<noradCross<<"xs="<<xs<<endl;
            //xs=xs*(d_E*d_omiga/num_evt);  //in unit of mub now
            xs=noradCross*(d_E*d_omiga/num_evt);  //in unit of mub now
            //cout<<"d_omiga="<<d_omiga<<"d_E="<<d_E<<endl;
            rate = xs * 1.0e-6 * 1e-24 * lumi;   //in unit of Hz
            //cout<<"norad rate="<<rate<<"xs="<<xs<<endl;
            
            //Jixie: add doRadiateBornXS to control this block from command line argument, sometimes we  
            //skip radiating the born xs to speed up
            if(rad_status>0 && doRadiateBornXS!=0){
                double RL_before_all=0;
                if(rad_status==1) RL_before_all=RL_before+(vz-vertex_z_min)/(vertex_z_max-vertex_z_min)*RL;
                if(rad_status==2) RL_before_all=RL_before+0.5*RL;
                //Jixie: since Eloss has been considered, RL_before and RL_after need to be set to zero
                RL_after=RL_before_all=0.0;

                // cout<<(vz-vertex_z_min)/(vertex_z_max-vertex_z_min) << " " <<RL_before_all<<endl;
                RadiativeCorrections rad;
                rad.SetTa(RL_after);
                rad.SetTb(RL_before_all);
                rad.SetCrossSection(noXS);
                double radCross=rad.Radiate();
                // cout<<"radCross="<<radCross<<endl;
                raddXSdEdOmega_mubGeVSr=radCross;
                radxs=radCross*(d_E*d_omiga/num_evt);  //in unit of mub now
                raterad = radxs * 1.0e-6 * 1e-24 * lumi;   //in unit of Hz
                // cout<<"rad rate="<<rate<<"radxs="<<radxs<<endl;
            }
            //calculate PVDIS asymmetries Abeam and AL
            //proton
            double A_beam_proton=calculate_proton_Abeam(unpol_pdf,  x,  Q2,  y);
            double A_L_proton=calculate_proton_AL(unpol_pdf,  pol_pdf, x,  Q2,  y);
            //neutron
            double A_beam_neutron=calculate_neutron_Abeam(unpol_pdf,  x,  Q2,  y);
            double A_L_neutron=calculate_neutron_AL(unpol_pdf,  pol_pdf, x,  Q2,  y);
            // do an average for the asymmetry according to Z and A
            Abeam= (Z*A_beam_proton + (A-Z)*A_beam_neutron)/A;
            AL= (Z*A_L_proton + (A-Z)*A_L_neutron)/A;


            T->Fill();
            count++;

            //output to lund file
            OUTPUT_lund << "1" << " \t " << x << " \t " << y  << " \t " << W << " \t " << Q2 << " \t " << rate << " \t " << raterad  << " \t " << 0  << " \t " << Abeam <<" \t " << AL << " \t " << Ei << " \t " << xi << " \t " << yi << " \t " << Wi << " \t " << Q2i << endl;

            //output to lund file (old format)
            //OUTPUT_lund << "1" << " \t " << Abeam  << " \t " << AL  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << rate << endl;

            OUTPUT_lund << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ep << " \t " << mass<<" \t " << vx  << " \t " << vy << " \t " << vz << endl;
        }  // save something physical for DIS events
        
        if(((i+1)%1000)==0 || i+1==num_evt) {
            cout<<i+1<<"/"<<num_evt<<" events thrown, "<<count<<" events saved. \r";
        }
    } //end for

// 	int counter_try=0,counter_good=0;
// 	for(int i=0; i<1e10; i++){
// 		if(counter_good>=num_evt) break;
// 		counter_try++;
//
// 		//uniform vx, vy, vz
// 		vx=rand.Uniform(vertex_x_min, vertex_x_max);
// 		vy=rand.Uniform(vertex_y_min, vertex_y_max);
// 		vz=rand.Uniform(vertex_z_min, vertex_z_max);
//
// 		//---------phase space
// 		//uniform Ep
// 		Ep=rand.Uniform(Ep_min, Ep_max);
// 		//uniform phi, theta
// 		phi=rand.Uniform(0,360);  //0 degree to 360 degree
// 		theta=rand.Uniform(theta_min, theta_max);
//
//
// 		double d_omiga=2*PI*(cos(theta_min*deg_to_rad) - cos(theta_max*deg_to_rad));
// 		double d_E=Ep_max-Ep_min;
//
// 		px=Ep*sin(theta*deg_to_rad)*cos(phi*deg_to_rad);
// 		py=Ep*sin(theta*deg_to_rad)*sin(phi*deg_to_rad);
// 		pz=Ep*cos(theta*deg_to_rad);
//
// 		//calculate kinematics
// 		double Nu=E-Ep;
// 		Q2=4.0*E*Ep*sin(theta*deg_to_rad/2.0)*sin(theta*deg_to_rad/2.0);
// 		W=sqrt(proton_mass*proton_mass + 2*proton_mass*Nu - Q2);
// 		x=Q2/2/proton_mass/Nu;
// 		y=Nu/E;
//
//
// 		if(x>=0 && x<=1){
// 			counter_good++;
//
// 			xs=calculate_fixed_target_xs( E,  Z,  A,  theta,  Ep,  unpol_pdf);   //theta unit in degree
// 			//xs in unit of mub/GeV-sr
// 			xs=xs*(d_E*d_omiga);  //in unit of mub now
//
// 			rate_pre = xs * 1.0e-6 * 1e-24 * lumi;   //in unit of Hz
//
// 			//calculate PVDIS asymmetries Abeam and AL
// 			//proton
// 			double A_beam_proton=calculate_proton_Abeam(unpol_pdf,  x,  Q2,  y);
// 			double A_L_proton=calculate_proton_AL(unpol_pdf,  pol_pdf, x,  Q2,  y);
// 			//neutron
// 			double A_beam_neutron=calculate_neutron_Abeam(unpol_pdf,  x,  Q2,  y);
// 			double A_L_neutron=calculate_neutron_AL(unpol_pdf,  pol_pdf, x,  Q2,  y);
// 			// do an average for the asymmetry according to Z and A
// 			Abeam= (Z*A_beam_proton + (A-Z)*A_beam_neutron)/A;
// 			AL= (Z*A_L_proton + (A-Z)*A_L_neutron)/A;
//
// 		T->Fill();
//
// 		}  // save something physical for DIS events
// 	} //end for
//
// TBranch *brate = T->Branch("rate", &rate, "data/D");       //in unit Hz
//
// // 			cout << counter_try << endl;
// 		for(int i=0; i<num_evt; i++){
// 		  T->GetEntry(i);
// // 		  cout << rate_pre << "  ";
// 		  rate=rate_pre/double(counter_try);
// // 		  cout << rate << "  ";
// 		  brate->Fill();
// 		//output to lund file
// 		OUTPUT_lund << "1" << " \t " << Abeam  << " \t " << AL  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << rate << endl;
// 		OUTPUT_lund << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ep << " \t " << mass<<" \t " << vx  << " \t " << vy << " \t " << vz << endl;
// 		}
// 		cout << endl;
//
// /*		for(int i=0; i<num_evt; i++){
// 		  T->GetEntry(i);
// 		  cout << rate_pre << "  ";
// 		  cout << rate << "  ";
// 		}
// 		cout << endl;	*/

    myfile->Write();
    myfile->Close();

    return 0;

}
