#ifndef EVGEN_INCLUSIVE_E_INPUT_PARAMETERS_HH
#define EVGEN_INCLUSIVE_E_INPUT_PARAMETERS_HH

// data structure for input parameters 

#include <cstdlib>
#include <string>

typedef struct inputParameters {
   std::string output_name;       // output file name  
   std::string pol_pdfset_name;   // polarized PDF model name
   std::string pol_pdfset_ID;     // polarized PDF model ID 
   std::string unpol_pdfset_name; // unpolarized PDF model name
   std::string unpol_pdfset_ID;   // unpolarized PDF model ID  
   double lumi;                   // [cm^-2 s^-1] per NUCLEUS, not NUCLEON 
   double E_beam;                 // [GeV] beam energy  
   double theta_min;              // [deg] scattering angle (min) 
   double theta_max;              // [deg] scattering angle (max) 
   double Ep_min;                 // [GeV] scattered electron energy (min) 
   double Ep_max;                 // [GeV] scattered electron energy (max) 
   double vx_min;                 // [cm]  vertex information.  depends on raster size  
   double vx_max;                 // [cm]  vertex information.  depends on raster size  
   double vy_min;                 // [cm]  vertex information.  depends on raster size 
   double vy_max;                 // [cm]  vertex information.  depends on raster size 
   double vz_min;                 // [cm]  vertex information. 
   double vz_max;                 // [cm]  vertex information. 
   int num_evt;                   // number of events  
   int tgt_Z;                     // target Z: number of protons       
   int tgt_A;                     // target A: atomic mass (= Z + N)
   int scale;                     // 0: no scale factor, 1:apply a scale factor on the cross sections
   int rad;                     // 0: Born cross section, 1:radiative cross section to get rate
   double RLb;                     //Total radiation length X0 before vertex 
   double RLa;                     //Total radiation length X0 after vertex 

   // constructor 
   inputParameters(): 
      output_name("NONE"),pol_pdfset_name("NONE"),pol_pdfset_ID("0000"),unpol_pdfset_name("NONE"),unpol_pdfset_ID("0000"),
      lumi(0),E_beam(0),theta_min(0),theta_max(0),Ep_min(0),Ep_max(0),
      vx_min(0),vx_max(0),vy_min(0),vy_max(0),vz_min(0),vz_max(0),
      num_evt(0),tgt_Z(0),tgt_A(0),scale(0),rad(0),RLb(0),RLa(0)
   {}

} inputParameters_t; 

#endif 
