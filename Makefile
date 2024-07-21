
LLVM_PATH= /opt/homebrew/Cellar/llvm/18.1.8/#Default path, please change here (versioning might be enough)
OPENBLAS_DIR = /opt/homebrew/Cellar/openblas/0.3.27/#Default path, please change here (versioning might be enough)
OPENMP_DIR = /opt/Homebrew/Cellar/libomp/18.1.6#Default path, please change here (versioning might be enough)
ARM_PER_LIB = /opt/arm/armpl_24.04_flang-new_clang_18#Default path, please change here (versioning might be enough)


CURRENT_DIR = $(shell pwd)
EXEC = tests
CXXFLAGS =  -Xclang -std=c++20 -fopenmp=libomp -O3 -flto -Wall
VARS = -DMATRIXMARKET -DMACOS
CXX = $(LLVM_PATH)bin/clang++
CC = $(LLVM_PATH)bin/clang


INCLUDES = -Iinclude/ \
	-I$(OPENBLAS_DIR)\
	-I$(OPENMP_DIR)/lib \
	-I$(ARM_PER_LIB)/include
	
LIBRSB =  $(shell echo $$HOME)/local/librsb/bin/
$(shell export PATH=PATH:$(LIBRSB))


SRC_PATH = src
LIB_DIRS = -L$(OPENBLAS_DIR)/lib  -L/usr/local/lib -L$(OPENMP_DIR)/lib -L$(ARM_PER_LIB)/lib
LIBS = -lopenblasp-r0.3.27  
#To activate ARMPL
#Uncomment line 30 and add -DARMPL in line 42, 43, 45
#LIBS = -lopenblasp-r0.3.27  -larmpl

SRCS = $(SRC_PATH)/main.cpp \
 $(SRC_PATH)/sparmatsymblk.cc \
 $(SRC_PATH)/MatSymBMtInstance.cpp \
 $(SRC_PATH)/matsym.cc
VERBOSE = -DVERBOSE=2
OBJS = $(SRCS:.cpp=.o)
OBJS := $(OBJS:.cc=.o)
TARGET = $(EXEC)
BLOCKSIZE = 3
macOS:
	$(MAKE) all VARS="-DMACOS -DCYTOSIM -DCYTOSIM_ORIGINAL -DCYTOSIM_NEW -DCYTOSIM_TEST -DBLOCKSIZE=$(BLOCKSIZE) $(VERBOSE)"
macOS-mmkt:
	$(MAKE) all VARS="-DMATRIXMARKET -DMACOS -DCYTOSIM -DCYTOSIM_ORIGINAL  -DCYTOSIM_TEST -DCYTOSIM_NEW -DBLOCKSIZE=$(BLOCKSIZE) $(VERBOSE)"
macOS-cytosimMat:
	$(MAKE) all VARS="-DCYTMAT -DMACOS -DCYTOSIM -DCYTOSIM_ORIGINAL -DCYTOSIM_NEW -DCYTOSIM_TEST -DBLOCKSIZE=$(BLOCKSIZE) $(VERBOSE)"
linux:
	$(MAKE) -f Makefile.linux VARS="-DCYTOSIM -DCYTOSIM_ORIGINAL -DCYTOSIM_NEW -DRSB -DBLOCKSIZE=$(BLOCKSIZE) $(VERBOSE)"
linux-mmkt:
	$(MAKE) -f Makefile.linux VARS="-DMATRIXMARKET -DCYTOSIM -DCYTOSIM_ORIGINAL -DCYTOSIM_NEW -DCYTOSIM_TEST -DRSB -DBLOCKSIZE=$(BLOCKSIZE) $(VERBOSE)"
linux-cytosimMat:
	$(MAKE) -f Makefile.linux VARS="-DCYTMAT -DCYTOSIM -DCYTOSIM_ORIGINAL -DCYTOSIM_NEW -DCYTOSIM_TEST -DRSB -DBLOCKSIZE=$(BLOCKSIZE) $(VERBOSE)"
linux-debug:
	$(MAKE) -f Makefile.linux VARS="-DMATRIXMARKET -DCYTOSIM_ORIGINAL -DBLOCKSIZE=$(BLOCKSIZE) $(VERBOSE)"
unvalid:
	$(MAKE) all VARS="-DCYTMAT -DMACOS -DMATRIXMARKET -DCYTOSIM -DCYTOSIM_ORIGINAL -DCYTOSIM_NEW -DCYTOSIM_TEST -DBLOCKSIZE=$(BLOCKSIZE) $(VERBOSE)"
all:
	$(MAKE) clean
	$(MAKE) compile
	$(MAKE) install
	
compile: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(VARS) $(LIB_DIRS)  -o $(TARGET) $(OBJS) $(LIBS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(VARS) -c $< -o $@

%.o : %.cc
	$(CXX) $(CCFLAGS) $(INCLUDES) $(VARS) -c $< -o $@

profile:
	$(MAKE) clean
	$(MAKE) all CXXFLAGS="$(CXXFLAGS) -pg -O3"

	
install:
	install_name_tool -add_rpath /usr/local/lib $(EXEC)
	install_name_tool -add_rpath $(OPENMP_DIR)/lib $(EXEC)
	install_name_tool -add_rpath $(ARM_PER_LIB)/lib $(EXEC)
clean:
	rm -f $(OBJS) $(TARGET)

