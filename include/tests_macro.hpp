#ifndef tests_macro_hpp
    #include <string>
    #include <iostream>
    #include "real.h"
    #ifndef VERBOSE
        #define VERBOSE 2 
    #endif 
    #define tests_macro_hpp
    #if VERBOSE <=1
        #define VLOG(x) ((void) 0)
        #define VLOGe(x) ((void) 0)
    #endif 
    #if VERBOSE >=2
        #define VLOG(x) std::cout <<" [INFO] "<<x << std::endl
        #define VLOGe(x) std::cerr <<" [ERROR] "<<x << std::endl
    #endif 

    template <typename T>
    inline std::string toStringPers(const T& value) {
        return std::to_string(value);
    }
    template <>
    inline std::string toStringPers(const std::string& value) {
        return value;
    }
    #define CONV(x) toStringPers(x)
    #define CLOCK std::chrono::high_resolution_clock::now();
    #define COUNT_T(start,stop) std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
    
  
    #define OUTPATH "res/compute.out"
    #define SEPARATOR " "
    #define FILE_HEADER "ARMPL RSB CYTOSIM_ORIGINAL CYTOSIM_NEW CYTOSIM_TEST\n";
    #define OPEN_OUTFILE std::ofstream outfile;\
            outfile.open(OUTPATH);\
            outfile<<FILE_HEADER;
    #define LOG_VALUE(start, stop) outfile<<COUNT_T(start, stop)<<SEPARATOR
    #define LOG_NULL outfile<<"-1 "
    #define CLOSE_OUTFILE outfile.close();

    #define DIFF_FILE_P "diffile.out"
    #define DIFF_FILE_HEADER "This file presents a detail of computation differences between algorithms\n\n"
    #if VERBOSE >=1
        #define DIFF_FILE std::ofstream diffile;\
                diffile.open(DIFF_FILE_P);\
                diffile<<DIFF_FILE_HEADER;
        
        #define LOG_VEC(Vector,Size,Name) diffile<<"\nValues computed for algorithm : "<<Name<<"\n";\
        for(int i=0; i<Size;i++){diffile<<Vector[i]<<" ";};

        #define EPSI 1
        inline int diff_count(real* Vector1, real* Vector2, int Size)
        {
            int nbDiff =0;
            for(int i=0; i<Size; i++)
            {
                if(abs(Vector1[i]-Vector2[i])>EPSI)
                {
                    nbDiff++;
                }
            }
            return nbDiff;
        }

        #define COUNT_DIF(Vector1,Vector2,Size) diff_count(Vector1, Vector2, Size)
        #define DIFF_FILE_C diffile.close();
        #define LOG_DIFF(Vector1, Vector2, Size,Name1,Name2) std::cout<<"Number of diff "<<Name1<<" - "<<Name2<<" : "<<CONV(COUNT_DIF(Vector1, Vector2, Size))<<std::endl;
    #endif
    #if VERBOSE == 0
        #define DIFF_FILE 
        #define LOG_DIFF(a,b,c,d,e)
        #define LOG_VEC(Vector,Size,Name)
        #define DIFF_FILE_C 
    #endif
    #endif