#!/usr/bin/env bash


omega_SM_ac.opt -scatter "e- e+ -> e- nuebar u dbar" -target:module cc20_ac_amplitude -target:whizard -target:parameter_module parameters_sm_ac > cc20_ac_amplitude.f90
gfortran -g -I$(whizard-config --prefix)/lib/mod/omega -I$(whizard-config --prefix)/lib/mod/whizard -I$(whizard-config --prefix)/lib/mod/models -fPIC -c cc20_ac_amplitude.f90 #-L$(whizard-config --prefix)/lib -lomega_core
g++ -g -Wl,-soname,libAmplitude.so -shared -o libcc20_ac_amplitude.so cc20_ac_amplitude.o -L$(whizard-config --prefix)/lib -lgfortran -lomega_core -lomega
