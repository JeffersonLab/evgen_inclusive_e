#ifndef ELOSS_PVDIS_H
#define ELOSS_PVDIS_H

#include "ELoss.h"
using namespace ELOSS;
using namespace std;

namespace PVDIS {

    static const int kNMatl = 5; //number of materials

    static int gTgtType = 1;  //1) LH2 target or 2) LD2 target
    static double gTgtEntrZPos = -10.0 *cm;  //the entrance face z location of the target

    static Material_t gMatl[kNMatl]; 
    ////////////////////////////////////////////////////////////////

    //construct the target geometry, in cm
    void ConstructTarget(int pTgtType, double pTgtEntrZPos);

    //calculate Mean energy loss for proton, in MeV and cm
    double CalELoss_proton(double E0, double vz, double *Ef);

    //calculate energy loss for electron, in MeV and cm
    double CalELoss_electron(double E0, double vz, int WithIntRad, double *Ef);
    double CalELoss_electron(double E0, double vz, int WithIntRad=0);


    //plot energy loss figure, in MeV
    void PlotELoss(double beamE, int nThrown);

}

#endif //ELOSS_PVDIS_H
