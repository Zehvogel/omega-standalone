#!/usr/bin/env bash


omega_SM.opt -scatter "e- e+ -> e- nuebar u dbar"  -target:parameter_module parameters_sm > cc20_amplitude.f90
gfortran -g -I$(whizard-config --prefix)/lib/mod/omega -fPIC -c parameters_SM.f90 cc20_amplitude.f90 #-L$(whizard-config --prefix)/lib -lomega_core
g++ -g -Wl,-soname,libAmplitude.so -shared -o libAmplitude.so parameters_SM.o cc20_amplitude.o -L$(whizard-config --prefix)/lib -lgfortran -lomega_core
