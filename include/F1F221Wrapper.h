// -*- C++ -*-

/* class F1F221Wrapper
 * Wrapper class to get f1, f2, and xs with P. Bosted fitting.
 * Unit of xs is ub/MeV-sr.
 * Valid for all W<5 GeV and all Q2<11 GeV2.
 * Use -fno-leading-underscore -fno-second-underscore when compiling F1F221.f
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   may 2014, Comments modified By Jixie Zhang

#ifndef F1F221WRAPPER_H
#define F1F221WRAPPER_H

class F1F221Wrapper {
public:
    F1F221Wrapper();
    ~F1F221Wrapper();

    void GetF1F2IN21(double Z, double A, double Q2, double W2, double &F1, double &F2);
    void GetF1F2QE21(double Z, double A, double Q2, double W2, double &F1, double &F2);
    double GetXS(double Z, double A, double Ei, double Ef, double theta);
};

#endif
