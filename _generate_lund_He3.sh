#!/bin/tcsh -f 
source /w/halla-scifs17exp/solid/tianye/solid_simulation/solid/evgen/evgen_inclusive_e-dflay_dev-4/setup
 /w/halla-scifs17exp/solid/tianye/solid_simulation/solid/evgen/evgen_inclusive_e-dflay_dev-4/solid_inclusive_e /w/halla-scifs17exp/solid/tianye/solid_simulation/solid/evgen/evgen_inclusive_e-dflay_dev-4/input/parameters_d2n.txt
  mv gen_d2n.root ./gen_d2n_fileN.root
  rm gen_d2n.lund

