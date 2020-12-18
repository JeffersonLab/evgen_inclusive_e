#!/bin/tcsh -f 
source /w/halla-scifs17exp/solid/tianye/solid_simulation/evgen_inclusive_e/setup
 /w/halla-scifs17exp/solid/tianye/solid_simulation/evgen_inclusive_e/solid_inclusive_e /w/halla-scifs17exp/solid/tianye/solid_simulation/evgen_inclusive_e/input/parameters_d2n.txt
  mv gen_d2n.root ./gen_d2n_fileN.root
  rm gen_d2n.lund

