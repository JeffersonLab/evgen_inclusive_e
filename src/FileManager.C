#include "FileManager.h"
#include <iterator>
using namespace std;
//_____________________________________________________________________________
FileManager::FileManager(){
   fVerbosity = 0;
}
//_____________________________________________________________________________
FileManager::~FileManager(){

}
//_____________________________________________________________________________
int FileManager::SplitString_whiteSpace(const std::string myStr,std::vector<std::string> &out){
   std::istringstream buffer(myStr);
   std::copy(std::istream_iterator<std::string>(buffer),
	 std::istream_iterator<std::string>(),
	 std::back_inserter(out));
   return 0;
}
//_____________________________________________________________________________
int FileManager::LoadInputData(const char *inpath,inputParameters_t &par){

   std::string aLine;
   std::vector<std::string> line;

   std::ifstream infile;
   infile.open(inpath);

   if( infile.fail() ){
      std::cout << "[FileManager::LoadInputData]: Cannot open the file: " << inpath << std::endl;
      return 1;
   }else{
      if(fVerbosity>0) std::cout << "[FileManager::ReadFile]: Opened the file: " << inpath << std::endl;
      while( !infile.eof() ){
	 std::getline(infile,aLine);
	 line.push_back(aLine);
      }
      line.pop_back();
      infile.close();
   }

   std::vector<std::string> data;
   // now parse the data 
   int NROW = line.size();
   for(int i=0;i<NROW;i++){
      // split to vector
      if(fVerbosity>0) std::cout << "[FileManager::LoadInputData]: Processing line: " << line[i] << std::endl; 
      SplitString_whiteSpace(line[i],data);
      // if line starts with a hash, skip (it's a comment) 
      if(data[0].compare("#")==0){
	 data.clear();
	 continue; 
      } 
      // now assign to the output container 
      if( data[0].compare("output_name")==0       ) par.output_name       = data[1]; 
      if( data[0].compare("pol_pdfset_name")==0   ) par.pol_pdfset_name   = data[1]; 
      if( data[0].compare("pol_pdfset_ID")==0     ) par.pol_pdfset_ID     = data[1]; 
      if( data[0].compare("unpol_pdfset_name")==0 ) par.unpol_pdfset_name = data[1]; 
      if( data[0].compare("unpol_pdfset_ID")==0   ) par.unpol_pdfset_ID   = data[1]; 
      if( data[0].compare("lumi")==0              ) par.lumi              = std::atof( data[1].c_str() ); 
      if( data[0].compare("E_beam")==0            ) par.E_beam            = std::atof( data[1].c_str() ); 
      if( data[0].compare("theta_min")==0         ) par.theta_min         = std::atof( data[1].c_str() ); 
      if( data[0].compare("theta_max")==0         ) par.theta_max         = std::atof( data[1].c_str() ); 
      if( data[0].compare("Ep_min")==0            ) par.Ep_min            = std::atof( data[1].c_str() ); 
      if( data[0].compare("Ep_max")==0            ) par.Ep_max            = std::atof( data[1].c_str() ); 
      if( data[0].compare("vx_min")==0            ) par.vx_min            = std::atof( data[1].c_str() ); 
      if( data[0].compare("vx_max")==0            ) par.vx_max            = std::atof( data[1].c_str() ); 
      if( data[0].compare("vy_min")==0            ) par.vy_min            = std::atof( data[1].c_str() ); 
      if( data[0].compare("vy_max")==0            ) par.vy_max            = std::atof( data[1].c_str() ); 
      if( data[0].compare("vz_min")==0            ) par.vz_min            = std::atof( data[1].c_str() ); 
      if( data[0].compare("vz_max")==0            ) par.vz_max            = std::atof( data[1].c_str() ); 
      if( data[0].compare("num_evt")==0           ) par.num_evt           = std::atoi( data[1].c_str() ); 
      if( data[0].compare("tgt_A")==0             ) par.tgt_A             = std::atoi( data[1].c_str() ); 
      if( data[0].compare("tgt_Z")==0             ) par.tgt_Z             = std::atoi( data[1].c_str() );
      // set up for next line of data 
      data.clear(); 
   }

   return 0;

}
