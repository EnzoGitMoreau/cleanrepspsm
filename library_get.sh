#!/bin/bash

wget https://sourceforge.net/projects/librsb/files/latest/download
mv download librsb.tar.gz
tar -xf librsb.tar.gz
cd librsb-1.3.0.2
./configure --enable-fortran-module-install --enable-matrix-types="double, double complex" \
	CC=gcc CXX=g++ FC=gfortran CFLAGS=-O3 CXXFLAGS=-O3 FCFLAGS=-O3 --prefix=$HOME/local