#ifndef EVGEN_INCLUSIVE_E_FILE_MANAGER_HH
#define EVGEN_INCLUSIVE_E_FILE_MANAGER_HH

// a file manager for evgen_inclusive_e parameters 

#include <cstdlib> 
#include <iostream>
#include <fstream>
#include <sstream>  
#include <vector>
#include <string>
#include <iterator>  

#include "InputParameters.h"

class FileManager {
   
   private:
      int fVerbosity; 

      int SplitString_whiteSpace(const std::string myStr,std::vector<std::string> &out); 

   public:
      FileManager();
      ~FileManager();

      void SetVerbosity(int v) { fVerbosity = v; } 

      int LoadInputData(const char *inpath,inputParameters_t &data);  

}; 

#endif 
