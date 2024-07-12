


//#define RSB 
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include "tests_macro.hpp"
//------------------------------Compilation checking------------------------------------------------
#ifndef BLOCKSIZE
    #define UNVALID
    #define ERROR_MESSAGE "Please specify a block_size"
#endif
#ifdef BLOCKSIZE
    #if BLOCKSIZE>=5
        #define UNVALID
        #define ERROR_MESSAGE "Block size must be 2,3 or 4"
    #endif 
    #if BLOCKSIZE <= 1
        #define UNVALID
        #define ERROR_MESSAGE "Block size must be 2,3 or 4"

    #endif 
#endif
#ifdef MATRIXMARKET
    #ifdef CYTMAT
        #ifndef UNVALID
        #define UNVALID
        #endif
        #ifndef ERROR_MESSAGE
        #define ERROR_MESSAGE "Unvalid compilation: Choose between MATRIXMARKET and CYTMAT, please refer to makefile"
        #endif
    #endif
#endif

//------------------------------End compilation checking------------------------------------------------
#ifdef MACOS
    #ifdef RSB
        #ifndef UNVALID
        #define UNVALID
        #endif 
        #ifndef ERROR_MESSAGE
        #define ERROR_MESSAGE "Unvalid compilation: no RSB possible on macOS"
        #endif 
    #endif 
#endif
#ifndef UNVALID
    #ifdef MATRIXMARKET
        #include "matrix_market_reader.hpp"
    #endif
    #ifdef CYTMAT
        #include "cytosim_matrix_reader.hpp"
    #endif
    #include "matsym.h"
    #include "sparmatsymblk.h"
    #include <random>
    #include "matrix.h"
    #include <time.h>
    #include "omp.h"
    #include <deque> 
    #include <tuple>
    #include "real.h"
    #include <chrono>
    #include <fstream>
    #include "custom_alg.h"

    #ifdef RSB
        #include <rsb.hpp>
    #endif
    #ifdef ARMPL
        #include "armpl.h"
        double* amd_matrix_vecmul(int size, int nTests, std::vector<std::pair<int, int>> pairs)
        {
            double* values = (double*) malloc(sizeof(double)*size*size);
            add_block_to_pos(values,pairs,size);
            
            
            armpl_spmat_t armpl_mat;
            armpl_int_t creation_flags = 0;
            armpl_status_t info = armpl_spmat_create_dense_d(&armpl_mat,ARMPL_COL_MAJOR,size, size, size, values,creation_flags);
            if (info!=ARMPL_STATUS_SUCCESS)
                printf("ERROR: armpl_spmat_create_csr_d returned %d\n", info);

            /* 3a. Supply any pertinent information that is known about the matrix */
            info = armpl_spmat_hint(armpl_mat, ARMPL_SPARSE_HINT_STRUCTURE,
                                ARMPL_SPARSE_STRUCTURE_UNSTRUCTURED);
            if (info!=ARMPL_STATUS_SUCCESS)
                printf("ERROR: armpl_spmat_hint returned %d\n", info);

            /* 3b. Supply any hints that are about the SpMV calculations
                to be performed */
            info = armpl_spmat_hint(armpl_mat, ARMPL_SPARSE_HINT_SPMV_OPERATION,
                                ARMPL_SPARSE_OPERATION_NOTRANS);
            if (info!=ARMPL_STATUS_SUCCESS)
                printf("ERROR: armpl_spmat_hint returned %d\n", info);

            info = armpl_spmat_hint(armpl_mat, ARMPL_SPARSE_HINT_SPMV_INVOCATIONS,
                                ARMPL_SPARSE_INVOCATIONS_MANY);
            if (info!=ARMPL_STATUS_SUCCESS)
                printf("ERROR: armpl_spmat_hint returned %d\n", info);

            /* 4. Call an optimization process that will learn from the hints you
                have previously supplied */
            info = armpl_spmv_optimize(armpl_mat);
            if (info!=ARMPL_STATUS_SUCCESS)
                printf("ERROR: armpl_spmv_optimize returned %d\n", info);
            double *x = (double *)malloc(size*sizeof(double));
            for (int i=0; i<size; i++) {
                    x[i] = i;
            }
            double *y = (double *)malloc(size*sizeof(double));
            for (int i=0; i<nTests; i++) {
                    info = armpl_spmv_exec_d(ARMPL_SPARSE_OPERATION_NOTRANS, 1.0,
                                    armpl_mat, x, 1.0, y);
                    if (info!=ARMPL_STATUS_SUCCESS)
                    printf("ERROR: armpl_spmv_exec_d returned %d\n", info);
            }

        
            armpl_spmat_destroy(armpl_mat);
            return y;
        }
    #endif 
    #ifdef RSB 
        void rsb_matrix_vecmul(double* Vec, double*res ,rsb::RsbMatrix<double>* mtx ,int nTests)
        {
            for (int i=0; i<nTests; i++)
            {
                mtx->spmv(RSB_TRANSPOSITION_N, 1.0, Vec, 1.0, res);
            }
        }
        rsb::RsbMatrix<double>* rsb_matrix_set(int size, std::vector<std::pair<int,int>> pairs)
        {
            rsb::RsbLib rsblib;
        const rsb_coo_idx_t nrA { size }, ncA { size };
        rsb::RsbMatrix<double>*mtx = new rsb::RsbMatrix<double>(nrA, ncA); // Declarations of IA,JA,VA are all accepted via <span>

        for (size_t k = 0; k < pairs.size(); k++) 
            { 
                int row = pairs[k].first;
                int col = pairs[k].second;
                for(int i=0; i<4; i++)
                {
                    for(int j=0; j<4;j++)
                    {
                        mtx->set_val((j/4.0), row*4 +i, col*4+j);
                    }
                }
            }
            mtx->close();
            return mtx;
        }
    #endif 


    int main(int argc, char* argv[])
    {
//------------------------------Base objects----------------------------------------------------------------
        using milli = std::chrono::milliseconds;
        auto start = CLOCK
        auto stop= CLOCK
        OPEN_OUTFILE;
        SparMatSymBlk testMatrix = SparMatSymBlk();
        
        int size = 0;
        int nb_threads;
        int nMatrix = 1;
        double block_percentage = 0.10;
//------------------------------Parsing user parameters-----------------------------------------------------
        #ifndef MATRIXMARKET
            #ifndef CYTMAT
                if (argc < 4)
                    {
                        VLOGe("Usage: tests nb_threads matSize");
                        return 1;
                    }
                else
                {
                    try
                    {
                        size = std::stoi(argv[2]);
                        nb_threads = std::stoi(argv[1]); 
                        nMatrix = std::stoi(argv[3]);
                    }
                    catch(std::exception e)
                    {
                        VLOGe("Usage: tests nb_threads matSize nMatrix [block_percentage]");
                        return 1;
                    }
                    if(size%(nb_threads*BLOCKSIZE)!= 0)
                    {
                        size = (int)(size / (BLOCKSIZE*nb_threads)) *BLOCKSIZE*nb_threads;
                    }
                    if(argc >=5)
                    {
                    try
                    {
                        block_percentage = std::stod(argv[4]);
                    }
                    catch(std::exception e)
                    {
                        VLOGe("block percentage not recognized as a double: usage x.f");
                        return 1;
                    }
                    }
                }
            #endif
        #endif 
        #ifdef MATRIXMARKET
            std::string inputPath;
            if(argc < 4)
            {
                VLOGe("Usage tests nbThreads nbMult matrixPath");
                return 1;
            }
            else
            {
                try
                {
                    nb_threads = std::stoi(argv[1]);
                    nMatrix = std::stoi(argv[2]);
                    inputPath = argv[3];
                }
                catch(const std::exception& e)
                {
                    VLOGe(e.what());
                    VLOGe("\n");
                }
                
            }
            #endif
            #ifdef CYTMAT
            std::string matPath;
            std::string vectorPath;
            if(argc < 5)
            {
                VLOGe("Usage tests nbThreads nbMult matrixPath vectorPath");
                return 1;
            }
            else
            {
                try
                {
                    nb_threads = std::stoi(argv[1]);
                    nMatrix = std::stoi(argv[2]);
                    matPath = argv[3];
                    vectorPath = argv[4];
                }
                catch(const std::exception& e)
                {
                    VLOGe(e.what());
                    VLOGe("\n");
                }
                
            }
        #endif
//------------------------------End parsing user parameters-------------------------------------------------
        omp_set_num_threads(nb_threads);
//------------------------------Generating random Matrix if matrix not given------------------------------------------------
        #ifndef MATRIXMARKET
            std::vector<std::pair<int, int>> pairs = select_random_points(size/4,(int) size*size/16 * block_percentage*block_percentage);
            testMatrix.resize(size);
            testMatrix.reset();
            #ifndef RSB
                add_block_to_pos_std(&testMatrix, pairs, size);
            #endif
            #ifdef RSB
            #ifndef CYTMAT
                rsb::RsbLib rsblib;
                rsb::RsbMatrix<double>* mtx;
                add_block_to_pos_rsb(&mtx, &testMatrix, pairs, size);
            #endif 
            #endif 
            VLOG("Constructed matrix of size "+CONV(size)<<" with "+CONV((size*block_percentage/BLOCKSIZE)*(size*block_percentage/BLOCKSIZE))+" blocks of size "+CONV(BLOCKSIZE));
            VLOG("Preparing to do"+CONV(nMatrix)+" multiplications");

        #endif
        
//------------------------------Reading given matrix------------------------------------------------
        #ifdef MATRIXMARKET
            VLOG("Asked to read "+CONV(inputPath)+" and to do "+CONV(nMatrix)+" mult with "+CONV(nb_threads)+" threads");
            #ifndef RSB 
                MatrixReader matrixReader(inputPath, &testMatrix,nb_threads);
            #endif 
            #ifdef RSB 
                rsb::RsbLib rsblib;
                rsb::RsbMatrix<double>* mtx;
                MatrixReader matrixReader(inputPath, &testMatrix,&mtx,nb_threads);
           
            #endif
            size = matrixReader.problemSize();
        #endif 
//------------------------------Reading given matrix + given vector------------------------------------------------
        #ifdef CYTMAT
            real* Vec;
            VLOG("Asked to read "+CONV(matPath)+" and "+CONV(vectorPath)+ " and to do "+CONV(nMatrix)+" mult with "+CONV(nb_threads)+" threads");
            #ifndef RSB 
                CytMatrixReader CytmatrixReader(matPath,vectorPath, &testMatrix,&Vec,nb_threads);
            #endif 
            #ifdef RSB 
            rsb::RsbLib rsblib;
            rsb::RsbMatrix<double>* mtx;
            
            CytMatrixReader CytmatrixReader(matPath, vectorPath, &testMatrix,&mtx, &Vec, nb_threads);
            #endif
            VLOG("Finished reading matrix and vector");
            size = CytmatrixReader.problemSize();
        #endif
//------------------------------Final inits------------------------------------------------
        
        real* Y_res = (real*)malloc(size * sizeof(real));//Y
        real* Y_true = (real*)malloc(size * sizeof(real));//Y_true to compare
        real* Y_rsb = (real*)malloc(size * sizeof(real));//Y_true to compare
        real* Y_dif = (real*)malloc(size * sizeof(real));//Y_diff that will store differences
        
        #ifndef CYTMAT
            real* Vec = (real*)malloc(size * sizeof(real));//Defining vector to do MX
        #endif
        #ifdef CYTOSIM_TEST 
            real* Y_test = (real*)malloc(size * sizeof(real));//Y_diff that will store differences
        #endif


        for(int i=0; i<size;i++)
        {
            Vec[i]=i;//Init
            Y_res[i] = 0;
            Y_true[i] = 0;
            Y_rsb[i] = 0;
        #ifdef CYTOSIM_TEST
            Y_test[i] = 0;
        #endif

        }
        testMatrix.prepareForMultiply(1);
        

//------------------------------Computing matrix multiplications for each algorithm------------------------------------------------
        VLOG("OpenMP is enabled with "+CONV(omp_get_max_threads())+" threads");
        #ifdef ARMPL
            real* y_arm = (real*)malloc(size * sizeof(real));//Y_true to compare
            VLOG("ARM PL is working");
            VLOG("Starting ARMPL matrix-vector multiplications");
            start = CLOCK
            y_arm = amd_matrix_vecmul(size, nMatrix, pairs);
            stop= STOP_T
            VLOG("ARMPL multiplications done in "+CONV(COUNT_T(start,stop))+" ms");
            outfile1 << COUNT_T(start,stop)<<" ";
        #endif 
        #ifndef ARMPL
        LOG_NULL;
        #endif
        #ifdef RSB
            VLOG("Starting libRsb matrix-vector multiplications");
            start = CLOCK
            rsb_matrix_vecmul(Vec, Y_rsb, mtx, nMatrix);
            stop= CLOCK
            LOG_VALUE(start, stop);
            VLOG("libRsb multiplications done in "+CONV(COUNT_T(start,stop))+" ms");
        #endif 
        #ifndef RSB
            LOG_NULL;
        #endif
        #ifdef CYTOSIM_ORIGINAL
            VLOG("Cytosim matrix-vector multiplications");
            start = CLOCK
            for(int i=0; i<nMatrix;i++)
            {
                testMatrix.vecMulAdd(Vec, Y_true);
            }
            stop= CLOCK
            LOG_VALUE(start, stop);
            VLOG("Cytosim multiplications done in "+CONV(COUNT_T(start,stop))+" ms");
        #endif 
        #ifdef CYTOSIM_TEST
            VLOG("Starting new-impl-test matrix-vector multiplications");
            start = CLOCK
            testMatrix.vecMulMtTest(nb_threads, Vec, Y_test,nMatrix);
            stop = CLOCK
            LOG_VALUE(start, stop);
            VLOG("new-impl-test multiplications done in "+CONV(COUNT_T(start,stop))+" ms");
        #endif
        #ifndef CYTOSIM_ORIGINAL
            LOG_NULL;
        #endif
        #ifdef CYTOSIM_NEW
            VLOG("Starting new-impl matrix-vector multiplications");
            start = CLOCK
                testMatrix.vecMulMt2(nb_threads, Vec, Y_res,nMatrix);
            stop = CLOCK
            LOG_VALUE(start, stop);
            VLOG("new-impl multiplications done in "+CONV(COUNT_T(start,stop))+" ms");
        #endif
        #ifndef CYTOSIM_NEW
            LOG_NULL;
        #endif
       
        CLOSE_OUTFILE;
        
//------------------------------Computing differences between algorithms results------------------------------------------------
        
        
        LOG_DIFF(Y_res, Y_true, size, "V1","Correct")
        LOG_DIFF(Y_test, Y_true, size, "V2","Correct")
        #ifdef RSB
            LOG_DIFF(Y_rsb, Y_true, size, "RSB","Correct")
        #endif
        DIFF_FILE
        LOG_VEC(Y_res, size, "Version 1")
        LOG_VEC(Y_true, size, "Cytosim")
        LOG_VEC(Y_test, size, "Version 2")
        #ifdef RSB 
            LOG_VEC(Y_rsb, size, "LibRSB")
        #endif
        
        DIFF_FILE_C
        
       
    }
#endif 
//------------------------------If compilation is wrong ( bad parameters )-----------------------------------------------
#ifdef UNVALID
    int main(int argc, char* argv[])
    {
        #ifdef ERROR_MESSAGE
            std::cout<<ERROR_MESSAGE;
        #endif
        #ifndef ERROR_MESSAGE
            std::cout<<"Unvalid compilation: no reason found";
        #endif
    }
#endif 