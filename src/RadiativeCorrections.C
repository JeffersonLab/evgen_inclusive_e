#include "RadiativeCorrections.h"
//________________________________________________________________________
RadiativeCorrections::RadiativeCorrections(){
   Init();
}
//________________________________________________________________________
RadiativeCorrections::~RadiativeCorrections(){

}
//________________________________________________________________________
void RadiativeCorrections::Init(){
   fDeltaE   = 0.01;           // in GeV
   fMT       = 0;
   fZ        = 0;
   fA        = 0;
   fb        = 0;
   fXi       = 0;
   fEta      = 0;
   fTa       = 0;
   fTb       = 0;
   fT        = 0;
   fThDeg    = 0;
   fEs       = 0;
   fEp       = 0;
   fR        = 0;
   fCFACT    = 0;
   fMT       = 0;
   fThreshold = RC::kPion; 
}
//_____________________________________________________________________________________________
double RadiativeCorrections::Radiate(){
   // set important variables 
   fZ     = fInclXS->GetZ();
   fA     = fInclXS->GetA();
   fEs    = fInclXS->GetEs();
   fEp    = fInclXS->GetEp();
   fThDeg = fInclXS->GetTh();
   fMT    = fA*proton_mass;      // set the target mass 
   fT     = fTa + fTb;

   if( (fTa==0)||(fTb==0) ){
      std::cout << "[RadiativeCorrections::Radiate]: WARNING! Radiation lengths are zero! " << std::endl;
   }

   CalculateEta();
   CalculateB();
   CalculateXi();
   CalculateR();
   CalculateCFACT();

   double BornXS = fInclXS->GetBornXS();
   double AnsEs  = CalculateEsIntegral();
   double AnsEp  = CalculateEpIntegral();
   double RadXS  = fCFACT*BornXS + AnsEs + AnsEp;

   return RadXS;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetPhi(double v){
   double phi = 1.0 - v + (3.0/4.0)*v*v;
   return phi;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetTr(double Q2){
   // General terms
   double M2 = electron_mass*electron_mass;
   // Individual terms
   double T1 = (1.0/fb)*(alpha/PI);
   double T2 = log(Q2/M2) - 1.0;
   // Put it all together 
   double Tr = T1*T2;
   return Tr;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetFTilde(double Q2){
   // General terms
   double M2     = electron_mass*electron_mass;
   double PI2    = PI*PI;
   double thr    = fThDeg*deg_to_rad;
   double COS    = cos(thr/2.0);
   double COS2   = COS*COS;
   double SPENCE = GetSpence(COS2);
   // Individual terms 
   double T1     = 1.0 + 0.5772*fb*fT;
   double T2     = (2.0*alpha/PI)*( (-14.0/9.0) + (13.0/12.0)*log(Q2/M2) );
   double T3     = (-1.0)*(alpha/(2.0*PI))*log( pow(fEs/fEp,2.0) );
   double T4     = (alpha/PI)*( (PI2/6.0) - SPENCE );
   // Put it all together
   double FTilde = T1+T2+T3+T4;
   return FTilde;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetEsMin(double Ep){
   double thr   = fThDeg*deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;

   double num=0,denom=0;
   if(fThreshold==RC::kElastic){
      num   = Ep;
      denom = 1.0 - (2.0*Ep/fMT)*SIN2; 
   }else if(fThreshold==RC::kPion){
      // this EXCLUDES the QE tail  
      num   = pion_mass*pion_mass + 2.*proton_mass*pion_mass + 2.*proton_mass*Ep;
      denom = 2.*proton_mass - (4.0*Ep)*SIN2;
   }

   double EsMin = num/denom;
   return EsMin;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetEpMax(double Es){
   double thr   = fThDeg*deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;

   double num=0,denom=0;
   if(fThreshold==RC::kElastic){
      num   = Es;
      denom = 1.0 + (2.0*Es/fMT)*SIN2; 
   }else if(fThreshold==RC::kPion){
      // this EXCLUDES the QE tail  
      num   = 2.*proton_mass*Es - 2.*proton_mass*pion_mass - pion_mass*pion_mass;
      denom = 2.*proton_mass + (4.0*Es)*SIN2;
   }

   double EpMax = num/denom;
   return EpMax;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::GetSpence(double x){
   // converted from radcor.f: 
   double num=0,denom=0,Index=0;
   double PI2 = PI*PI;
   double ans = (PI2/6.0) - log(x)*log(1.-x);

   for(int i=0;i<50;i++){
      Index   = (double)i + 1.0;
      num     = pow(x,i+1.0);
      denom   = pow(Index,2.);
      if(denom>0){
	 ans -= num/denom;
      }else{
	 ans -= 0;
      }
   }

   return ans;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateEta(){
   double Z23   = pow(fZ,-2.0/3.0);
   double Z13   = pow(fZ,-1.0/3.0);
   double num   = log(1440.0*Z23);
   double denom = log(183.0*Z13);
   fEta         = num/denom;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateB(){
   double Z13 = pow(fZ,-1.0/3.0);
   double T1  = 1.0;
   double T2  = (1.0/9.0)*( (fZ+1.0)/(fZ+fEta) );
   double T3  = 1.0/log(183.0*Z13);
   fb         = (4.0/3.0)*(T1 + T2*T3);
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateXi(){
   double Z13 = pow(fZ,-1.0/3.0);
   double T1  = PI*electron_mass/(2.0*alpha);
   double T2  = fT/( (fZ+fEta)*log(183.0*Z13) );
   fXi        = T1*T2;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateR(){
   double thr   = fThDeg*deg_to_rad;
   double SIN   = sin(thr/2.0);
   double SIN2  = SIN*SIN;
   double num   = fMT + 2.0*fEs*SIN2;
   double denom = fMT - 2.0*fEp*SIN2;
   fR           = num/denom;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::CalculateCFACT(){
   // General terms 
   double Q2     = Kinematics::GetQ2(fEs,fEp,fThDeg);
   double Tr     = GetTr(Q2);
   double FTilde = GetFTilde(Q2);
   // First term
   double Term1  = fR*fDeltaE/fEs;
   double Exp1   = fb*(fTb+Tr);
   double T1     = pow(Term1,Exp1);
   // Second term 
   double Term2  = fDeltaE/fEp;
   double Exp2   = fb*(fTa+Tr);
   double T2     = pow(Term2,Exp2);
   // Third term
   double num    = fXi/fDeltaE;
   double denom  = 1.0 - fb*(fTa+fTb+2.0*Tr);
   double T3     = 1.0 - num/denom;
   // Put it all together 
   fCFACT       = FTilde*T1*T2*T3;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::EsIntegrand(const double EsPrime){
   // general terms
   double Q2       = Kinematics::GetQ2(EsPrime,fEp,fThDeg);
   double FTilde   = GetFTilde(Q2);
   double Tr       = GetTr(Q2);
   double dEs      = fEs-EsPrime;
   double v        = dEs/fEs;
   double phi      = GetPhi(v);
   fInclXS->SetEs(EsPrime);
   fInclXS->SetEp(fEp);
   fInclXS->SetTh(fThDeg);
   double Sig      = fInclXS->GetBornXS();
   if(Sig!=Sig){
      std::cout << "[RadiativeCorrections::EsIntegrand]: Invalid cross section! " << std::endl;
      exit(1);
   }
   double SigTilde = FTilde*Sig;
   // first term 
   double Term1    = dEs/(fEp*fR);
   double Exp1     = fb*(fTa+Tr);
   double T1       = pow(Term1,Exp1);
   // second term
   double Term2    = dEs/fEs;
   double Exp2     = fb*(fTb+Tr);
   double T2       = pow(Term2,Exp2);
   // third term 
   double T3       = fb*( ((fTb+Tr)/dEs)*phi + fXi/(2.0*pow(dEs,2.0)) );
   // put it all together
   double FES      = T1*T2*T3*SigTilde;

   return FES;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::EpIntegrand(const double EpPrime){
   // general terms 
   double Q2       = Kinematics::GetQ2(fEs,EpPrime,fThDeg);
   double Tr       = GetTr(Q2);
   double FTilde   = GetFTilde(Q2);
   double dEp      = EpPrime-fEp;
   double v        = dEp/EpPrime;
   double phi      = GetPhi(v);
   fInclXS->SetEs(fEs);
   fInclXS->SetEp(EpPrime);
   fInclXS->SetTh(fThDeg);
   double Sig      = fInclXS->GetBornXS();
   if(Sig!=Sig){
      std::cout << "[RadiativeCorrections::EpIntegrand]: Invalid cross section! " << std::endl;
      exit(1);
   }
   double SigTilde = FTilde*Sig;
   // first term 
   double Term1    = dEp/(EpPrime);
   double Exp1     = fb*(fTa+Tr);
   double T1       = pow(Term1,Exp1);
   // second term
   double Term2    = (dEp*fR)/fEs;
   double Exp2     = fb*(fTb+Tr);
   double T2       = pow(Term2,Exp2);
   // third term 
   double T3       = fb*( ((fTa+Tr)/dEp)*phi + fXi/(2.0*pow(dEp,2.0)) );
   // put it all together
   double FEP      = T1*T2*T3*SigTilde;

   return FEP;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::CalculateEsIntegral(){
   int depth      = 10;
   double epsilon = 1e-10;
   double min     = GetEsMin(fEp);
   double max     = fEs - fR*fDeltaE;
   double AnsEs   = Integrate(&RadiativeCorrections::EsIntegrand,min,max,epsilon,depth);
   return AnsEs;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::CalculateEpIntegral(){
   int depth      = 10;
   double epsilon = 1e-10;
   double min     = fEp + fDeltaE;
   double max     = GetEpMax(fEs);
   double AnsEp   = Integrate(&RadiativeCorrections::EpIntegrand,min,max,epsilon,depth);
   return AnsEp;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::Integrate(double (RadiativeCorrections::*f)(const double),double A,double B,double epsilon,int Depth){
   // Adaptive Simpson's Rule
   double C   = (A + B)/2.0;
   double H   = B - A;
   double fa  = (this->*f)(A);
   double fb  = (this->*f)(B);
   double fc  = (this->*f)(C);
   double S   = (H/6.0)*(fa + 4.0*fc + fb);
   double ans = AdaptiveSimpsonAux(f,A,B,epsilon,S,fa,fb,fc,Depth);
   return ans;
}
//_____________________________________________________________________________________________
double RadiativeCorrections::AdaptiveSimpsonAux(double (RadiativeCorrections::*f)(const double),
      double A,double B,double epsilon,
      double S,double fa,double fb,double fc,int bottom){
   // Recursive auxiliary function for AdaptiveSimpson() function
   double C      = (A + B)/2.0;
   double H      = B - A;
   double D      = (A + C)/2.0;
   double E      = (C + B)/2.0;
   double fd     = (this->*f)(D);
   double fe     = (this->*f)(E);
   double Sleft  = (H/12.0)*(fa + 4.0*fd + fc);
   double Sright = (H/12.0)*(fc + 4.0*fe + fb);
   double S2     = Sleft + Sright;
   if( (bottom <= 0) || (fabs(S2 - S) <= 15.0*epsilon) ){
      return S2 + (S2 - S)/15;
   }
   double arg = AdaptiveSimpsonAux(f,A,C,epsilon/2.0,Sleft, fa,fc,fd,bottom-1) +
      AdaptiveSimpsonAux(f,C,B,epsilon/2.0,Sright,fc,fb,fe,bottom-1);
   return arg;
}
//_____________________________________________________________________________________________
void RadiativeCorrections::Print(){
   std::cout << "------------------------------------"              << std::endl;
   std::cout << "Radiative correction quantities: "                 << std::endl;
   std::cout << "DeltaE = " << std::fixed      << std::setprecision(4) << fDeltaE << " [GeV]"  << std::endl;
   std::cout << "Constants for given thicknesses: "                 << std::endl;
   std::cout << "Tb     = " << std::scientific << std::setprecision(4) << fTb     << " [#X0]"  << std::endl;
   std::cout << "Ta     = " << std::scientific << std::setprecision(4) << fTa     << " [#X0]"  << std::endl;
   std::cout << "eta    = " << std::scientific << std::setprecision(4) << fEta    << " [-]"    << std::endl;
   std::cout << "b      = " << std::fixed      << std::setprecision(4) << fb      << " [-]"    << std::endl;
   std::cout << "xi     = " << std::scientific << std::setprecision(4) << fXi     << " [GeV]"  << std::endl;
   std::cout << "Values that change for each (Es,Ep): "   << std::endl;
   std::cout << "Es         = " << std::fixed      << std::setprecision(4) << fEs     << " [GeV]"   << std::endl;
   std::cout << "Ep         = " << std::fixed      << std::setprecision(4) << fEp     << " [GeV]"   << std::endl;
   std::cout << "R          = " << std::fixed      << std::setprecision(4) << fR      << " [-]"     << std::endl;
   std::cout << "CFACT      = " << std::scientific << std::setprecision(4) << fCFACT  << " [-]"     << std::endl;
}
