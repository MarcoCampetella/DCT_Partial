#!/bin/sh

gfortran main.f90 util.f90 -fopenmp -O3 -funroll-loops -o ~/bin/dctSel.x
