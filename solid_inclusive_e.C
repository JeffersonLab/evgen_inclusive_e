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
  
char input_gen_file[50]; 
	    strcpy(input_gen_file,argv[1]);
    std::string str;
    std::string substring;
    str = input_gen_file;
    substring=str.substr(str.find("_")+1,3);
    std::cout<<"name"<< substring <<std::endl;
  //std::string outname(22,input_gen_file);
  //  std::cout << "input_file:        " << outname             << std::endl;
    // test the radiative correction class
    std::cout << "Creating RadiativeCorrection object..." << std::endl; 
    RadiativeCorrections *myRC = new RadiativeCorrections(); 
    std::cout << "--> Done!" << std::endl;
    delete myRC; 
    std::cout << "--> Deleted RC object" << std::endl;

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
    std::cout << "RLb:            " <<  par.RLb                  << std::endl;
    std::cout << "RLa:            " <<  par.RLa                  << std::endl;

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
	const int rad_status=par.rad;    //0 Born cross sections; 1:Radiative cross section to get rate  
	const double RL_before=par.RLb;   //radiation length before vertex 
	const double RL_after=par.RLa;   //radiation length after vertex 
	string pol_pdfset_name=par.pol_pdfset_name;    //pol. pdfset name
	string unpol_pdfset_name=par.unpol_pdfset_name;   // unpol. pdfset name
// 	string unpol_pdfset_name="CT14lo";   // unpol. pdfset name	
// 	string unpol_pdfset_name="cteq66";   // unpol. pdfset name		

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
// 	string unpol_pdfset_name="CT14lo";   // unpol. pdfset name	
// 	string unpol_pdfset_name="cteq66";   // unpol. pdfset name		
	TString name_rootfile_output=par.output_name;   // name of the rootfile to save output data
	
	//output lund file
	ofstream OUTPUT_lund;
	OUTPUT_lund.open(Form("gen_%s.lund",substring.data()));
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
        // radiative effects before vertex 
	cout<<"Tb is loaded......using "<<RL_before<<endl;
	cout<<"Ta is loaded......using "<<RL_after<<endl;
	

	//print some information
	vector<int> pol_pids=pol_pdf->flavors();
	for(int i=0; i<pol_pids.size();i++){
		 cout<<"# YX  : (pol pids) "<<i<<"	"<<pol_pids[i]<<endl;
	}
	
	vector<int> unpol_pids=unpol_pdf->flavors();
	for(int i=0; i<unpol_pids.size();i++){
		 cout<<"# YX  : (unpol pids) "<<i<<"	"<<unpol_pids[i]<<endl;
	}
	
	//##################################################################################
	//
	//
	//
	//               Define outputs
	//
	//##################################################################################               
	
	//TTree to save
	double Abeam=0, AL=0, x=0, y=0, W=0, Q2=0, rate=0,rate_pre=0;
	int charge=-1, particle_id=11;
	double px=0, py=0, pz=0;
	double Ep=0;
	double mass=0.511/1000.0;  //GeV
	double vx=0, vy=0, vz=0;
	double xs=0; 
	double dXSdEdOmega_mubGeVSr=0; // differential cross section 
	double radxs=0; 
	double raddXSdEdOmega_mubGeVSr=0; // differential cross section 
	double theta=0;
	double phi=0;
        double Nu=0;
        double factor=0;
	TFile *myfile=new TFile(name_rootfile_output,"RECREATE");
	TTree *T=new TTree("T","T");
	T->Branch("Abeam", &Abeam, "data/D");
	T->Branch("AL", &AL, "data/D");
	T->Branch("x", &x, "data/D");
	T->Branch("y", &y, "data/D");
	T->Branch("W", &W, "data/D");
	T->Branch("Q2", &Q2, "data/D");
	T->Branch("rate_pre", &rate_pre, "data/D");       //before normalized
	T->Branch("charge",&charge,"data/I");
	T->Branch("particle_id",&particle_id,"data/I");
	T->Branch("px",&px, "data/D");
	T->Branch("py",&py, "data/D");
	T->Branch("pz",&pz, "data/D");
	T->Branch("Ep",&Ep, "data/D");
	T->Branch("mass",&mass, "data/D");
	T->Branch("vx",&vx, "data/D");
	T->Branch("vy",&vy, "data/D");
	T->Branch("vz",&vz, "data/D");
	T->Branch("xs",&xs, "data/D");
	T->Branch("radxs",&radxs, "data/D");
	T->Branch("dXSdEdOmega_mubGeVSr",&dXSdEdOmega_mubGeVSr,"data/D");
	T->Branch("raddXSdEdOmega_mubGeVSr",&raddXSdEdOmega_mubGeVSr,"data/D");
	T->Branch("theta",&theta, "data/D");
	T->Branch("phi",&phi, "data/D");
	T->Branch("nu" ,&Nu, "data/D");

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

		//calculate kinematics
		Nu=E-Ep;
		Q2=4.0*E*Ep*sin(theta*deg_to_rad/2.0)*sin(theta*deg_to_rad/2.0);
		W=sqrt(proton_mass*proton_mass + 2*proton_mass*Nu - Q2);
		x=Q2/2/proton_mass/Nu;
		y=Nu/E;

		
		if(x>=0 && x<=1){
                         
			xs=calculate_fixed_target_xs( E,  Z,  A,  theta,  Ep,  unpol_pdf);   //theta unit in degree
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
                        NoRdXSDEdOmega.SetEs(E);
                        NoRdXSDEdOmega.SetEp(Ep);
                        NoRdXSDEdOmega.SetTh(theta);
                        NoRdXSDEdOmega.Setpdf(unpol_pdf);
                        NoRdXSDEdOmega.SetScale(factor);
                        // cout<<"dXSdEdOmega="<<dXSdEdOmega_mubGeVSr<<endl;
                         noXS=&NoRdXSDEdOmega;
                         double noradCross=noXS->GetBornXS();
			dXSdEdOmega_mubGeVSr = noradCross; 
                        // cout<<"noradCross="<<noradCross<<"xs="<<xs<<endl;
                         RadiativeCorrections rad;
                         rad.SetTa(RL_after);
                         rad.SetTb(RL_before); 
                         rad.SetCrossSection(noXS); 
                         double radCross=rad.Radiate();
                        // cout<<"radCross="<<radCross<<endl;
                         raddXSdEdOmega_mubGeVSr=radCross; 
			//xs=xs*(d_E*d_omiga/num_evt);  //in unit of mub now
			xs=noradCross*(d_E*d_omiga/num_evt);  //in unit of mub now
			radxs=radCross*(d_E*d_omiga/num_evt);  //in unit of mub now
                        //cout<<"d_omiga="<<d_omiga<<"d_E="<<d_E<<endl;
                        if(rad_status==0){	
			rate = xs * 1.0e-6 * 1e-24 * lumi;   //in unit of Hz
                         //cout<<"norad rate="<<rate<<"xs="<<xs<<endl;
                        }else{
			rate = radxs * 1.0e-6 * 1e-24 * lumi;   //in unit of Hz
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

		//output to lund file
		OUT << "1" << " \t " << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << rate << " \t " << 0  << " \t " << 0  << " \t "  << Abeam <<" \t " << AL << endl;
		
		//output to lund file (old format)
// 		OUTPUT_lund << "1" << " \t " << Abeam  << " \t " << AL  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << rate << endl;

		OUTPUT_lund << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ep << " \t " << mass<<" \t " << vx  << " \t " << vy << " \t " << vz << endl;
		}  // save something physical for DIS events
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
