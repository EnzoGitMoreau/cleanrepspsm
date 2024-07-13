
<h1>
Symmetric Sparse Matrix (block oriented) matrix vector multiplication algorithm
</h1>

## How to install

## Macos

1- Install Homebrew
2- install llvm's clang by with homebrew: brew install llvm
3- install openblas with homebrew : brew install openblas 
4- report installation directories in the makefile
5 (optional)- install ARMPL library to compare with ARM algorithms


## Linux

1- Download & install LibRSB
by executing:

 1. wget https://sourceforge.net/projects/librsb/files/latest/download
 2. mv download librsb.tar.gz
tar -xf librsb.tar.gz
cd librsb-(VERSION) *Enter your version here*
./configure --enable-fortran-module-install --enable-matrix-types="double, double complex" \
	CC=gcc CXX=g++ FC=gfortran CFLAGS=-O3 CXXFLAGS=-O3 FCFLAGS=-O3 --prefix=$HOME/local
	
 3. Add librsb to PATH by typing : 
 export PATH=\$PATH:$HOME/local/librsb/bin/
 Before each compilation

<h1>How to compile</h1>
<h2>macOS</h2>
<p>Requires ARMPL library for comparaison, install ARMPL on ARM website then compile with make macOS
<br>
Requires openBlas for matrix specifications, please install openBlas on openBlas website and indicate the headers directory in INCLUDE, and the lib directory in LIB.
<br>
To use a matrixMarket matrix in the .mtx format, please compile with make macOS-mmtk
<br>
To use a matrix in the .mtx format and a vector in .mtx format, please compile with make macOS-cytosimMat
</p>
<h2>Linux</h2>
<p>Requires libRSB for comparaison, install libRsb using libRsb repository.
Requires openBlas for matrix specifications, please install openBlas on openBlas website and indicate the headers directory in INCLUDE, and the lib directory in LIB
<h1>Usage</h1>
<h2>Normal mode</h2>
<p>tests [nbThreads] [matrixSize] [numberOfMultiplications] [blockDensity]<br>
Program will generate a random SparseSymmetric block matrix based on your specifications.<br>
Algorithms will do numberOfMultiplications Matrix-Vector multiplications and duration times will be printed<br>
To use a matrix in the .mtx format and a vector in .mtx format, please compile with make linux-mmkt<br>
To use a matrix in the .mtx format and a vector in .mtx format, please compile with make linux-cytosimMat<br>
</p>
<h2>MatrixMarket mode</h2>
<p>tests [nbThreads] [numberOfMultiplications] [pathToMatrix]<br>
Algorithms will do numberOfMultiplications Matrix-Vector multiplications and duration times will be printed
</p>

<h2>cytosimMarket mode</h2>
<p>tests [nbThreads] [numberOfMultiplications] [pathToMatrix] [pathToVector]<br>
Algorithms will do numberOfMultiplications Matrix-Vector multiplications and duration times will be printed
</p>

## Verbosity options

To manage verbosity options, refer to the VERBOSITY line in the Makefile.

 - -DVERBOSITY = 2, full logging of computation process, writing of matrix multiplications results and differences between algorithms and computation times.
 - -DVERBOSITY = 1, writing of matrix multiplications results and differences between algorithms and computation times.
- -DVERBOSITY = 0, writing of computation times.

Config file for output is include/tests_macro.hpp
OUT_PATH specifies the output for computation times.
 *Default is res/compute.out*
DIFFILE_P specifies the output for differences between algorithms. 
*Default is diffile.out*

## Important files

- in src: main.cpp generates tests.
- in include: 
	- MatSymBMtInstance.hpp implements the new algorithm for multi-threaded : SPSM Vec-Mul
	- matrix_market_reader.hpp and cytosim_matrix_reader.hpp implements parsing of matrix files.
	- custom_alg.hpp gives functions for random matrix generations 

