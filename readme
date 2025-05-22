# Introduction 
--------------------
This is an inclusive electron generator with different model and fit

F1F2 2021 fit
The W2<24 GeV2 and all Q2<30 GeV2---Eric Christy's fit (2021) significantly improved the fits to proton (W2>7 GeV2),deuteron, 4He, 12C, 27Al, and 56Fe (covering the full 12 GeV kinematics). For 3He, all the physics components in the fit scale very well with A but the Fermi momentum parameters has not been set to a reasonable values (plan to update the fit soon). 

F1F2 2009 fit + PDF
The process includes QE + resonance + DIS for nucleus target with Z protons and A atoms.
The W<3 GeV region uses Peter Bosted fit (https://userweb.jlab.org/~bosted/fits.html)
the W>3 GeV region uses world PDF sets, the LHAPDF6 interface is used
refer to http://hallaweb.jlab.org/12GeV/SoLID/download/sim/talk/Inclusive_electron_generator.pdf

Please send questions Zhiwen Zhao (zwzhao@jlab.org), Ye Tian (tainye@jlab.org)

# how to run it 
--------------------
refer to how to run pre-compiled version in the container at jlab ifarm 
https://github.com/JeffersonLab/solid_release/blob/master/howto_evgen.md

general running instruction as follows, don't use this
> source setup (for the container only)
>./evgen_inclusive_e <input_file> [doIncomingEloss] [doRadiateBornXS]

a root file and a txt file in lund format will be produced in the current directory

see inputfile examples under input dir
the main parameters are as follows:
* num_evt: number of trial events thrown evenly in costheta,phi,P space. only physical events in 0<x_bj<1 will be kept. thus output events are less than this
* Fit_model: 21 with 2021 fit, 9 with 2009 fit + pdf
* rad: 0 Born cross section only; 1:unpolarized radiative cross section with vertex z smearing; 2:unpolarized radiative cross section without vertex z smearing
* scale: 0 for without, 1 for with additional correction (0.906-0.00699*E) for F1F2_09 fit only, according to David Flay, it's need to compare to data
* ThreshType, Integration threshold set up for radiation corrections 0 for kElastic, 1 for kPion, and 2 for kQuasi-free
the optional command line options are
* doIncomingEloss=1 will turn on Eloss for beam before vertex, default is 0
* doRadiateBornXS=0 will turn rad off by ignoring input file, default is 1

# how to analyze output 
--------------------
setup env at ifarm
module use /group/halla/modulefiles
module load root/6.32.08

some output at /work/halla/solid/evgen/evgen_inclusive_e/

* plot eAll distribution with and without SoLID acceptance
> root 'analysis_eAll.C("eDIS_rootfile","acceptance_rootfile",beam_energy)'
> root 'analysis_eAll.C("/work/halla/solid/evgen/evgen_inclusive_e/commit0a256eb_20250521/eAll_solid_SIDIS_He3_11GeV_theta5-35deg_fit21_rad0_1e6.root","/group/solid/solid_github/JeffersonLab/solid_gemc/analysis/acceptance/result_SIDIS_He3/201701_24GeV/acceptance_solid_SIDIS_He3_electron_1e7_201701_output_final.root",11)'
> root 'analysis_eAll.C("/work/halla/solid/evgen/evgen_inclusive_e/commit0a256eb_20250521/eAll_solid_SIDIS_NH3_11GeV_theta3-90deg_fit21_rad0_5e6.root","/group/solid/solid_github/JeffersonLab/solid_gemc/analysis/acceptance/result_SIDIS_NH3/202012_24GeV/acceptance_solid_SIDIS_NH3_24GeV_electron_1e7_202012_output_final.root",11)'

 config                                   rate gen(khz)      rate accepted from both FA and LA(khz)
 solid_SIDIS_He3_8.8GeV_theta5-35deg      1165             185        
 solid_SIDIS_He3_11GeV_theta5-35deg       792              118
 solid_SIDIS_He3_22GeV_theta5-35deg       226              25
 solid_SIDIS_NH3_8.8GeV_theta3-90deg      565              2.7        
 solid_SIDIS_NH3_11GeV_theta3-90deg       379              1.75
 solid_SIDIS_NH3_22GeV_theta3-90deg       114              0.44

* compare this generator, nicked named "eAll", and "eDIS" generator (https://github.com/JeffersonLab/evgen_inclusive)
> root 'compare_eDIS_eAll.C("eDIS_rootfile","eAll_rootfile")'

* compare rad and norad
> root analysis_rad_norad_2D_plot.C

# log of major changes 
--------------------
2025/5/21 commit 2f12a68
let scale only apply to model 9 and ignore it for fit model 21, by Zhiwen Zhao

2023/9/8 commit d4898a1 
fix normalization bug by Ye Tian 
modified evgen_inclusive_e.cxx to throw events based on random costheta rather than theta (because in the calculate_fixed_target_xs function, the cross sections are calcualted based on costheta). 
**the corresponding change of the PVDIS rate could increase the uncertainties of the submited proposal that used the previous eAll generator (10-20 %).

2022/05/25 commit 6fc41a0
add ThreshType option for radiative correction to speed up calculation, by Ye Tian

2021/04/06 commit 4bb274c
finished adding F1F2 fit 2021 (model 21) in additional to old F1F2 fit 2009 (model 9) by Ye Tian and add beam energy after external rad in output by Jixie Zhang

2021/03/16 commit 2972a2e
modified by Jixie Zhang 
1) Add energy loss calculation routines: PVDIS::CalELoss_electron(beam, vz)
2) After Ep generated flatly in range (0,beam), call PVDIS::CalELoss_electron(beam, vz)  to get Ei, which is the energy of the incoming beam electron. If (Ei < Ep), call  PVDIS::CalELoss_electron () to get another Ei, repeat this till Ei > Ep,
3) change all beam energy to Ei when calculate Q2, nu,y and XS
4) change Ta and Tb to zero before calling David's subroutine to radiated born XS. Ta set to zero since energy loss have been taken care of by my Eloss routines, Tb is set to zero since Geant4 will take care of it. David's code will do the internal radiation correction to the born XS.
5) Add Ei into output root tree, and store it in UserVar008 in the lund file
6) Add usage to main(), print out usage if no argument is given 
7) Print out how many events have been saved when each 100 events thrown
8) Modify CMakelist, also add my own makefile to compile the program.  Both ways can compile the program

2021/01/11 commit 11d7751
adding unpolarized radiative correction and input file by David Flay and Ye Tian, output adjustment by Zhiwen Zhao

2020/04/24 commit 4b22e1b
convert to github by Zhiwen Zhao from https://jlabsvn.jlab.org/svnroot/solid/evgen/solid_inclusive_e

2017/02/13 initial version
written and tested with F1F2 fit 2009 and LHAPDF6.1.6 by Yuxiang Zhao (yxzhao@jlab.org)

# how to install
--------------------
refer to how to install it using the container at jlab ifarm 
https://github.com/JeffersonLab/solid_release/blob/master/howto_evgen.md

It requires ROOT5.34.36 and LHAPDF6.1.6 with pdfsets of NNPDFpol11_100 and CT14nlo at least
(see how to install LHAPDF6.1.6 below)

The current setup file has default ROOT and LHAPDF6 path on ifarm
You need to put your ROOT and LHAPDF6 path if compiling on other machine

>source setup
>cmake .
>make

a temporary compiling tool, don't use this
make -f Makefile_jx

install LHAPDF6:
refer to https://lhapdf.hepforge.org/install.html
for 6.1.6
> make sure boost packages are installed by checking "rpm -qa |grep boost"
> mkdir -p LHAPDF/source
> cd LHAPDF/source
> wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.1.6.tar.gz
> tar xf LHAPDF-6.1.6.tar.gz
> cd LHAPDF-6.1.6
> ./configure --prefix=$PWD/../../LHAPDF-6.1.6
> make
> make install
> cd ../../LHAPDF-6.1.6/share/LHAPDF/
> wget http://www.hepforge.org/archive/lhapdf/pdfsets/v6.backup/6.1.6/NNPDFpol11_100.tar.gz
> tar zxf NNPDFpol11_100.tar.gz
> wget http://www.hepforge.org/archive/lhapdf/pdfsets/v6.backup/6.1.6/CT14nlo.tar.gz
> tar zxf CT14nlo.tar.gz
> wget http://www.hepforge.org/archive/lhapdf/pdfsets/v6.backup/6.1.6/CT14lo.tar.gz
> tar zxf CT14lo.tar.gz
> wget http://www.hepforge.org/archive/lhapdf/pdfsets/v6.backup/6.1.6/cteq66.tar.gz
> tar zxf cteq66.tar.gz


# code details 
--------------------
The main code "evgen_inclusive_e.C" is segmented into several regions:
1. user inputs for the simulation -> define the simulation parameters
2. load polpodf and unpolpdfs for W>3 GeV region ->user doesn't need to deal with it normally
3. Define outputs, lund output and root tree definition  -> user can modify it according to the needs
4. throw events and calculate rate in src/fixed_target_xs.cxx, kinematics etc. for each event using the functions provided by the package   -->user can modify it according to the needs

All the functions are shown in the include directory, here is also a list of functions
that can be used by users:

//function from Peter Bosted and Eric Christy's electron fit
GetF1F2IN21(Z, A, Q2, W2, F1, F2)
F1F2IN09(int Z, int IA, double qsq, double wsq, double &F1, double &F2, double &Rc);

//function from LHAPDF
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
//neutron structure functions
double calculate_neutron_g1gz(PDF* pol_pdf, double x, double Q2);
double calculate_neutron_g5gz(PDF* pol_pdf, double x, double Q2);
double calculate_neutron_F2g(PDF* unpol_pdf, double x, double Q2);
double calculate_neutron_F1g(PDF* unpol_pdf, double x, double Q2);
double calculate_neutron_F1gz(PDF* unpol_pdf, double x, double Q2);
double calculate_neutron_F3gz(PDF* unpol_pdf, double x, double Q2);
//PVDIS asymmetry using long. pol. electron + unpolarized proton
double calculate_proton_Abeam(PDF* unpol_pdf, double x, double Q2, double y);
double calculate_neutron_Abeam(PDF* unpol_pdf, double x, double Q2, double y);
double calculate_fixed_target_xs(double E, int Z, int A, double theta, double Ep, PDF* unpol_pdf);
//PVDIS asymmetry using unpol electron + longitudinally polarized proton
double calculate_proton_AL(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_proton_AL_g1gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_proton_AL_g5gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_neutron_AL(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_neutron_AL_g1gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
double calculate_neutron_AL_g5gz(PDF* unpol_pdf, PDF* pol_pdf, double x, double Q2, double y);
