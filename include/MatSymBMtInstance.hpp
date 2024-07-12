//
//  MatSymBMtInstance.hpp
//  MatrixCalculation
//
//  Created by Moreau Enzo on 17/04/2024.
//

#ifndef MatSymBMtInstance_hpp
#define MatSymBMtInstance_hpp
#ifndef BLOCKSIZE 
#define FORMAT

#endif


#include <barrier>
#include <fstream>
#include <thread>
#include <stdio.h>
#include <mutex>
#include <chrono>
#include "sparmatsymblk.h"
#include "real.h"



class MatSymBMtInstance final
{
private:
    int phase_number;
    std::vector<int> valid_phase_index;
    int* phase_tab;
    int nbThreads;
    std::mutex _mutex;
    int big_mat_size;
    real* Y1;
    real* Y2;
    const real* X;
    SparMatSymBlk* matrix;
    int thNb =0;
    int threadBlockSize;
    real*** workingPhases;
    real** newPhases;
    int** newIndexes;
    int*** workingIndexes;
   

    int** work_lengths;
public:
    MatSymBMtInstance(SparMatSymBlk* matrix, int nbThreadsE)
    {
        nbThreads = nbThreadsE;
        this->matrix = matrix;
        big_mat_size = matrix->size();
        real* big_Y = (real*) malloc(big_mat_size*2*sizeof(real));
        Y1 = big_Y;
        Y2 = big_Y + big_mat_size;
        for(int i =0; i<big_mat_size;i++)
        {
            Y1[i] = 0;
            Y2[i] = 0;
        }
    }
    void generateBlocks()
    {
      
        phase_tab = (int*) malloc(sizeof(int)*nbThreads);
        phase_tab[0] = 0;
        if(nbThreads %2 == 0)
        {
            phase_number = nbThreads /2 +1;
            if(nbThreads == 4)
            {
                phase_tab[0] = 0;
                phase_tab[1] = 1;
                phase_tab[2] = 2;
                phase_tab[3] = 1;

            }
            else
            {
            for(int i =1; i<=nbThreads/2;i++)
            {
                phase_tab[i] = i;
                
            }
            for(int i=1 ; i<nbThreads/2;i++)
            {
                phase_tab[nbThreads/2 + i] = nbThreads/2 - i;
              
            }
            }
        }
        else
        {
          
            phase_number  = (nbThreads +1) / 2;
            if(nbThreads == 3)
            {
                phase_tab[0] = 0;
                phase_tab[1] = 1;
                phase_tab[2] = 1;
            }
            if(nbThreads == 5)
            {
                phase_tab[0] = 0;
                phase_tab[1] = 1;
                phase_tab[2] = 2;
                phase_tab[3] = 2;
                phase_tab[4] = 1;

            }
            if(nbThreads == 7)
            {
                phase_tab[0] = 0;
                phase_tab[1] = 1;
                phase_tab[2] = 2;
                phase_tab[3] = 3;
                phase_tab[4] = 3;
                phase_tab[5] = 2;
                phase_tab[6] = 1;

            }
             if(nbThreads == 9)
            {
                phase_tab[0] = 0;
                phase_tab[1] = 1;
                phase_tab[2] = 2;
                phase_tab[3] = 3;
                phase_tab[4] = 4;
                phase_tab[5] = 4;
                phase_tab[6] = 3;
                phase_tab[7] = 2;
                phase_tab[8] = 1;

            }
             if(nbThreads == 11)
            {
                phase_tab[0] = 0;
                phase_tab[1] = 1;
                phase_tab[2] = 2;
                phase_tab[3] = 3;
                phase_tab[4] = 4;
                phase_tab[5] = 5;
                phase_tab[6] = 5;
                phase_tab[7] = 4;
                phase_tab[8] = 3;
                phase_tab[9] = 2;
                phase_tab[10] = 1;

            }
        }
        
        workingPhases = (real***) malloc(sizeof(real**)*phase_number);
        workingIndexes = (int***)malloc(sizeof(int**)*phase_number);
        
        QueueBlock*** phasesTemp = (QueueBlock***) malloc(sizeof(QueueBlock**) * phase_number);
        threadBlockSize = (int) matrix->size() / (BLOCKSIZE * nbThreads);
        int* phase_counter = (int*) malloc(sizeof(int)*phase_number);
        work_lengths = (int**) malloc(sizeof(int*)*phase_number);
        for(int i=0; i<phase_number;i++)
        {
            work_lengths[i] = (int*) malloc(sizeof(int)*nbThreads);
        }
        for(int i =0; i<phase_number; i++)
        {
            workingPhases[i] = (real**) malloc(sizeof(real*) * nbThreads);
            workingIndexes[i] = (int**) malloc(sizeof(int*) * nbThreads);
            phasesTemp[i] = (QueueBlock**) malloc(sizeof(QueueBlock*)*nbThreads);
            for(int j=0; j<nbThreads; j++)
            {
                phasesTemp[i][j] = new QueueBlock();
            }
        }
        int** map_threads = (int**) malloc(nbThreads * sizeof(int*));
        for(int i = 0;i<nbThreads; i++)
        {
            map_threads[i] = (int*) malloc(nbThreads * sizeof(int));
        }

        for(int k =0; k<=nbThreads/2; k++)
            {
            for(int i =0; i<nbThreads;i++)
                {
                map_threads[(i+k)%nbThreads][i] = i;
                map_threads[i][(i+k)%nbThreads] = i;     
                }
            }
        #ifdef DEBUG
            std::cout<<"[DEBUG]Printing wanted matrix\n";
                        for(int i =0; i<nbThreads; i++)
                    {
                        for(int j=0; j<nbThreads; j++)
                        {
                            std::cout<<map_threads[i][j]<<" ";
                        }
                        std::cout<<"\n";
                    }
        #endif 
        int i =0;
        while(true)//Fetching blocks, assigned them to a phase each
        {
            //std::cout<<"here\n";
            int indX, indY;
            SparMatSymBlk::Column* col = &matrix->column_[i];
            if(col->nbb_>0)
            {
                indY = matrix->colidx_[i]/threadBlockSize;
                int indice_col = matrix->colidx_[i];
                
                for(int j=0; j<col->nbb_; j++)
                {
                    int indice_ligne = col->inx_[j];
                    indX =col->inx_[j]/threadBlockSize;
                    int phase =phase_tab[indX-indY];
                    Block data;
                    //std::cout<<"there\n";
                    #if BLOCKSIZE ==4
                    Matrix44 matrix= col->blk_[j];
                    #endif
                    #if BLOCKSIZE == 3
                    Matrix33 matrix = col->blk_[j];
                    #endif
                    data.a0 = matrix.val[0x0];
                    data.a1 = matrix.val[0x1];
                    data.a2 = matrix.val[0x2];
                    data.a3 = matrix.val[0x3];
                    data.a4 = matrix.val[0x4];
                    data.a5 = matrix.val[0x5];
                    data.a6 = matrix.val[0x6];
                    data.a7 = matrix.val[0x7];
                    data.a8 = matrix.val[0x8];
                    #if BLOCKSIZE ==4
                    data.a9 = matrix.val[0x9];
                    data.aA = matrix.val[0xA];
                    data.aB = matrix.val[0xB];
                    data.aC = matrix.val[0xC];
                    data.aD = matrix.val[0xD];
                    data.aE = matrix.val[0xE];
                    data.aF = matrix.val[0xF];
                    #endif 
                    data.index_x = indice_ligne * BLOCKSIZE;
                    data.index_y = indice_col * BLOCKSIZE;
                    //std::cout<<"la\n";
                    //std::cout<<"Matrix Size:"<<big_mat_size<<"\n";
                    //std::cout<<"threadblocksize:"<<threadBlockSize<<"\n";
                    //std::cout<<"ix"<<data.index_x<<"iy:"<<data.index_y<<"\n";
                    //std::cout<<"les donnees sont: "<<indX<<","<<indY<<"phase:"<<phase<<"\n";
                    //std::cout<<"et puis:"<<map_threads[indX][indY]<<"\n";
                    
                    phasesTemp[phase][map_threads[indX][indY]]->push_front(data);
                    //std::cout<<"j'ai reussi celui ci indX ne peut pas valoir 4 \n";
                    phase_counter[phase]++;
                    //std::cout<<"et pq pas la\n";
                }
            }
            //std::cout<<"carrement ici\n";
            i = matrix->colidx_[i+1];
            if(i>=matrix->rsize_)
            {
                
                break;
            }
            
        }
        //std::cout<<"Work lengths:\n";
        for(int i=0; i<phase_number; i++)
        {
            //std::cout<<"\nPhase: "<<i<<"\n";
            for(int j=0; j<nbThreads;j++)
            {
                
                
                int t_work_length = (int)phasesTemp[i][j]->size();
                //std::cout<<t_work_length<<" ";
                work_lengths[i][j] = t_work_length;
                
                workingPhases[i][j] = (real*) malloc(sizeof(real)*t_work_length*BLOCKSIZE*BLOCKSIZE);
                workingIndexes[i][j] = (int*) malloc(sizeof(int)*t_work_length*2);
                
                
                for(int m=0; m<t_work_length; m++)
                {
                    
                   
                    Block first= phasesTemp[i][j]->front();
                    phasesTemp[i][j]->pop_front();
                    
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m] = first.a0;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+1] = first.a1;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+2] = first.a2;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+3] = first.a3;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+4] = first.a4;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+5] = first.a5;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+6] = first.a6;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+7] = first.a7;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+8] = first.a8;
                    #if BLOCKSIZE == 4
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+9] = first.a9;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+10] = first.aA;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+11] = first.aB;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+12] = first.aC;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+13] = first.aD;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+14] = first.aE;
                    workingPhases[i][j][BLOCKSIZE*BLOCKSIZE*m+15] = first.aF;
                    #endif 
                    workingIndexes[i][j][2*m] = first.index_x;
                    workingIndexes[i][j][2*m+1] = first.index_y;
              
                    
                }
                
                
            }
        }

    
        for(int i =0; i<phase_number; i++)
        {
           
            int count = 0;
            
            for(int thNb = 0; thNb<nbThreads; thNb++)
            {
                count += work_lengths[i][thNb];
            }
            if(count>0)
            {
                valid_phase_index.emplace_back(i);
            }
        }
        
    }   
    int thread_number()
    {
        thNb++;
        return thNb-1;
    }
    void work2(std::barrier<>& barrier,std::barrier<>& barrier2, int n_work)
    {
        _mutex.lock();
        int thread_nb = thread_number();
        _mutex.unlock();
        
        for(int m = 0; m<n_work;m++)
        {
            workThread2(barrier2, thread_nb);
            barrier.arrive_and_wait();
        }
    } 
    void work4(std::barrier<>& barrier,std::barrier<>& barrier2, int n_work)
    {
       
        _mutex.lock();
        int thread_nb = thread_number();
        _mutex.unlock();
        
        for(int m = 0; m<n_work;m++)
        {
            
            workThread4(barrier2, thread_nb);
            
        }
    } 
    void workThread2(std::barrier<>& barrier2, int thNb)
        {
        

        int active_work_count;
        int j = 0;
        
        real* work;
        int work_nb;
        
        const real* X1;
        const real* X2;
        real* Y1_;
        real* Y2_;
        bool swap = false;
        while(j<valid_phase_index.size())
        {
            
            barrier2.arrive_and_wait();

            int k  = valid_phase_index[j];
            
            swap = thNb >= (nbThreads - k);
           
            work =  workingPhases[k][thNb]; 
            
            work_nb = work_lengths[k][thNb];
            int* index = workingIndexes[k][thNb];
            
            
            for(int m = 0; m<work_nb; m++)
            {
                int ix = index[m*2+0];
                int iy = index[m*2+1];
                real a0 = work[BLOCKSIZE*BLOCKSIZE*m+0];
                real a1 = work[BLOCKSIZE*BLOCKSIZE*m+1];
                real a2 = work[BLOCKSIZE*BLOCKSIZE*m+2];
                real a3 = work[BLOCKSIZE*BLOCKSIZE*m+3];
                real a4 = work[BLOCKSIZE*BLOCKSIZE*m+4];
                real a5 = work[BLOCKSIZE*BLOCKSIZE*m+5];
                real a6 = work[BLOCKSIZE*BLOCKSIZE*m+6];
                real a7 = work[BLOCKSIZE*BLOCKSIZE*m+7];
                real a8 = work[BLOCKSIZE*BLOCKSIZE*m+8];
                #if BLOCKSIZE ==3 
                real a9 = work[BLOCKSIZE*BLOCKSIZE*m+9];
                real aA = work[BLOCKSIZE*BLOCKSIZE*m+10];
                real aB = work[BLOCKSIZE*BLOCKSIZE*m+11];
                real aC = work[BLOCKSIZE*BLOCKSIZE*m+12];
                real aD = work[BLOCKSIZE*BLOCKSIZE*m+13];
                real aE = work[BLOCKSIZE*BLOCKSIZE*m+14];
                real aF = work[BLOCKSIZE*BLOCKSIZE*m+15];
                #endif 
                if(swap)
                {
                    X1 = X  + iy; //no swap
                    X2 = X  + ix;
                    Y1_= Y1 + ix;
                    Y2_= Y2 + iy;
                    
                }
                else
                {
                    X1 = X  + iy;
                    X2 = X  + ix;
                    Y1_= Y2 + ix;
                    Y2_= Y1 + iy;
                    
                    
                }
                
                if(ix-iy !=0)
                {  
                    #if BLOCKSIZE==4
                        real y10 = 0;
                        real y11 = 0;
                        real y12=0;
                        real y13=0;
                        real y20 = 0;
                        real y21 = 0;
                        real y22 = 0;
                        real y23=0;
                        real r10 = X1[0];
                        real r11 = X1[1];
                        real r12 = X1[2];
                        real r13 = X1[3];
                        real r20 = X2[0];
                        real r21 = X2[1];
                        real r22 = X2[2];
                        real r23 = X2[3];
                        {real & c = a0;
                            
                            y10 += c*r10;
                            y20 += c*r20;
                        }
                        {
                            real & c = a1;
                            y11 += c * r10;
                            y20 += c * r21;
                        }
                        {
                            real & c = a4;
                            y10 += c * r11;
                            y21 += c * r20;
                        }
                        {
                            real & c = a2;
                            y12 += c * r10;
                            y20 += c * r22;
                        }
                        {
                            real& c = a8;
                            y10 += c * r12;
                            y22 += c* r20;
                        }
                        {
                            real & c = a3;
                            y13 += c * r10;
                            y20 += c * r23;
                        }
                        {
                            real & c  = aC;
                            y10 += c *r13;
                            y23 += c * r20;
                        }
                        {
                            real & c = a5;
                            y11+= c*r11;
                            y21 += c *r21;
                        }
                        {
                            real & c=  a6;
                            y12 += c*r11;
                            y21 += c*r22;
                        }
                        {
                            real & c = a9;
                            y11 += c *r12;
                            y22 += c *r21;
                        }
                        {
                            real &c  = a7;
                            y13 += c *r11;
                            y21 += c * r23;
                        }
                        {
                            real & c  = aD;
                            y11 +=c*r13;
                            y23 += c *r21;
                        }
                        {
                            real & c  = aA;
                            y12 += c * r12;
                            y22 += c* r22;
                        }
                        {
                            real & c = aE;
                            y12 += c * r13;
                            y23 += c * r22;
                        }
                        {
                            real & c = aB;
                            y13 += c * r12;
                            y22 += c *r23;
                        }
                        {
                            real & c = aF;
                            y13 += c * r13;
                            y23 += c *r23;
                        }
                        
                        //_mutex.lock();
                        Y1_[0] += y10;
                        Y1_[1] += y11;
                        Y1_[2] += y12;
                        Y1_[3] += y13;
                        Y2_[0] += y20;
                        Y2_[1] += y21;
                        Y2_[2] += y22;
                        Y2_[3] += y23;
                        //_mutex.unlock();
                    #endif 
                    #if BLOCKSIZE==3
                    real y10 = 0;
                    real y11 = 0;
                    real y12=  0;
                    real y20 = 0;
                    real y21 = 0;
                    real y22 = 0;
                    real r10 = X1[0];
                    real r11 = X1[1];
                    real r12 = X1[2];
                    real r20 = X2[0];
                    real r21 = X2[1];
                    real r22 = X2[2];

                        {real & c = a0;
                            
                            y10 += c*r10;
                            y20 += c*r20;
                        }
                        {
                            real & c = a1;
                            y11 += c * r10;
                            y20 += c * r21;
                        }
                        {
                            real & c = a3;
                            y10 += c * r11;
                            y21 += c * r20;
                        }
                        {
                            real & c = a2;
                            y12 += c * r10;
                            y20 += c * r22;
                        }
                        {
                            real& c = a6;
                            y10 += c * r12;
                            y22 += c* r20;
                        }
                        {
                            real & c = a4;
                            y11 += c * r11;
                            y21 += c * r21;
                        }
                        {
                            real & c = a7;
                            y11+= c*r12;
                            y22 += c *r21;
                        }
                        {
                            real & c = a5;
                            y12 += c *r11;
                            y21 += c *r22;
                        }
                        {
                            real &c  = a8;
                            y12 += c *r12;
                            y22 += c * r22;
                        }
                    Y1_[0] += y10;
                    Y1_[1] += y11;
                    Y1_[2] += y12;
                    Y2_[0] += y20;
                    Y2_[1] += y21;
                    Y2_[2] += y22;
                    
                        
                    #endif 
                    #if BLOCKSIZE==2
                    #endif   
                }
                else
                {
                    
                    #if BLOCKSIZE == 4
                        Y1_[0] += a0 * X1[0] + a4 * X1[1] + a8 * X1[2] + aC * X1[3];
                        Y1_[1] += a1 * X1[0] + a5 * X1[1] + a9 * X1[2] + aD * X1[3];
                        Y1_[2] += a2 * X1[0] + a6 * X1[1] + aA * X1[2] + aE * X1[3];
                        Y1_[3] += a3 * X1[0] + a7 * X1[1] + aB * X1[2] + aF * X1[3];
                    #endif 
                    #if BLOCKSIZE == 3 
                        Y1_[0] += a0 * X1[0] + a3 * X1[1] + a6 * X1[2];
                        Y1_[1] += a1 * X1[0] + a4 * X1[1] + a7 * X1[2];
                        Y1_[2] += a2 * X1[0] + a5 * X1[1] + a8 * X1[2];
                     
                    #endif 
                    
                }
            }
            
           
            j++;
            barrier2.arrive_and_wait();
            
        }
        }
    void workThread4(std::barrier<>& barrier2, int thNb)
        {
            int active_work_count;
            int j = 0;
            real* work;
            int work_nb;
            const real* X1;
            const real* X2;
            real* Y1_;
            int ix;
            int iy;
            int ix_p;
            int iy_p;
            real* Y_base;
            real* Y_swap;
            real* Y2_;
            bool swap;
            while(j<valid_phase_index.size())
            {
              
                real acc_10 = 0;
                real acc_11 = 0;
                real acc_12 = 0;
                
                real acc_00 = 0;
                real acc_01 = 0;
                real acc_02 = 0;
                #if BLOCKSIZE == 4
                    real acc_03 = 0;
                    real acc_13 = 0;
                #endif 
                int k  = valid_phase_index[j];
                work =  workingPhases[k][thNb];
                swap = thNb >= (nbThreads - k);
                work_nb = work_lengths[k][thNb];
                int* index = workingIndexes[k][thNb];
                X1 = X;
                X2 = X;
                if(swap)
                {
                    Y_base = Y2;
                    Y_swap = Y1;
                 
                }
                else
                {   
                    Y_base = Y1;
                    Y_swap = Y2;
              
                }
                if(work_nb >0)
                {
                    
                    ix_p = index[0];
                    iy_p = index[1];
                    ix = ix_p;
                    iy = iy_p;
                    real a0 = work[0];
                    real a1 = work[1];
                    real a2 = work[2];
                    real a3 = work[3];
                    real a4 = work[4];
                    real a5 = work[5];
                    real a6 = work[6];
                    real a7 = work[7];
                    real a8 = work[8];
                    #if BLOCKSIZE == 4
                    real a9 = work[9];
                    real aA = work[10];
                    real aB = work[11];
                    real aC = work[12];
                    real aD = work[13];
                    real aE = work[14];
                    real aF = work[15];
                    #endif 
                    X1= X+iy;
                    X2= X+ix;
                    Y1_= Y_base + ix;
                    Y2_= Y_swap + iy;
                    //std::cout<<"thread : "<<thNb<<"starting with block"<<"("<<ix<<","<<iy<<")\n";
                    if(ix-iy !=0)
                        {
                            
                            #if BLOCKSIZE==4
                            real y10 = 0;
                            real y11 = 0;
                            real y12=0;
                            real y13=0;
                            real y20 = 0;
                            real y21 = 0;
                            real y22 = 0;
                            real y23=0;
                            real r10 = X1[0];
                            real r11 = X1[1];
                            real r12 = X1[2];
                            real r13 = X1[3];
                            real r20 = X2[0];
                            real r21 = X2[1];
                            real r22 = X2[2];
                            real r23 = X2[3];
                            {real & c = a0;
                                
                                y10 += c*r10;
                                y20 += c*r20;
                            }
                            {
                                real & c = a1;
                                y11 += c * r10;
                                y20 += c * r21;
                            }
                            {
                                real & c = a4;
                                y10 += c * r11;
                                y21 += c * r20;
                            }
                            {
                                real & c = a2;
                                y12 += c * r10;
                                y20 += c * r22;
                            }
                            {
                                real& c = a8;
                                y10 += c * r12;
                                y22 += c* r20;
                            }
                            {
                                real & c = a3;
                                y13 += c * r10;
                                y20 += c * r23;
                            }
                            {
                                real & c  = aC;
                                y10 += c *r13;
                                y23 += c * r20;
                            }
                            {
                                real & c = a5;
                                y11+= c*r11;
                                y21 += c *r21;
                            }
                            {
                                real & c=  a6;
                                y12 += c*r11;
                                y21 += c*r22;
                            }
                            {
                                real & c = a9;
                                y11 += c *r12;
                                y22 += c *r21;
                            }
                            {
                                real &c  = a7;
                                y13 += c *r11;
                                y21 += c * r23;
                            }
                            {
                                real & c  = aD;
                                y11 +=c*r13;
                                y23 += c *r21;
                            }
                            {
                                real & c  = aA;
                                y12 += c * r12;
                                y22 += c* r22;
                            }
                            {
                                real & c = aE;
                                y12 += c * r13;
                                y23 += c * r22;
                            }
                            {
                                real & c = aB;
                                y13 += c * r12;
                                y22 += c *r23;
                            }
                            {
                                real & c = aF;
                                y13 += c * r13;
                                y23 += c *r23;
                            }
                            
                            
                            Y1_[0] += y10;
                            Y1_[1] += y11;
                            Y1_[2] += y12;
                            Y1_[3] += y13;
                            Y2_[0] += y20;
                            Y2_[1] += y21;
                            Y2_[2] += y22;
                            Y2_[3] += y23;
                            
                            #endif 
                            #if BLOCKSIZE==3
                                real y10 = 0;
                                real y11 = 0;
                                real y12=  0;
                                real y20 = 0;
                                real y21 = 0;
                                real y22 = 0;
                                real r10 = X1[0];
                                real r11 = X1[1];
                                real r12 = X1[2];
                                real r20 = X2[0];
                                real r21 = X2[1];
                                real r22 = X2[2];

                                    {real & c = a0;
                                        
                                        y10 += c*r10;
                                        y20 += c*r20;
                                    }
                                    {
                                        real & c = a1;
                                        y11 += c * r10;
                                        y20 += c * r21;
                                    }
                                    {
                                        real & c = a3;
                                        y10 += c * r11;
                                        y21 += c * r20;
                                    }
                                    {
                                        real & c = a2;
                                        y12 += c * r10;
                                        y20 += c * r22;
                                    }
                                    {
                                        real& c = a6;
                                        y10 += c * r12;
                                        y22 += c* r20;
                                    }
                                    {
                                        real & c = a4;
                                        y11 += c * r11;
                                        y21 += c * r21;
                                    }
                                    {
                                        real & c = a7;
                                        y11+= c*r12;
                                        y22 += c *r21;
                                    }
                                    {
                                        real & c = a5;
                                        y12 += c *r11;
                                        y21 += c *r22;
                                    }
                                    {
                                        real &c  = a8;
                                        y12 += c *r12;
                                        y22 += c * r22;
                                    }
                                Y1_[0] += y10;
                                Y1_[1] += y11;
                                Y1_[2] += y12;
                                Y2_[0] += y20;
                                Y2_[1] += y21;
                                Y2_[2] += y22;
                                    
                            #endif 
                            #if BLOCKSIZE==2
                            //to_add
                            #endif 
                            

                        
                        }
                    else
                        {
                          
                            #if BLOCKSIZE ==4
                            acc_10 += a0 * X1[0] + a4 * X1[1] + a8 * X1[2] + aC * X1[3];
                            acc_11 += a1 * X1[0] + a5 * X1[1] + a9 * X1[2] + aD * X1[3];
                            acc_12 += a2 * X1[0] + a6 * X1[1] + aA * X1[2] + aE * X1[3];
                            acc_13 += a3 * X1[0] + a7 * X1[1] + aB * X1[2] + aF * X1[3];
                            #endif
                            #if BLOCKSIZE == 3
                            acc_10 += a0 * X1[0] + a3 * X1[1] + a6 * X1[2];
                            acc_11 += a1 * X1[0] + a4 * X1[1] + a7 * X1[2];
                            acc_12 += a2 * X1[0] + a5 * X1[1] + a8 * X1[2];
                            #endif 
                        
                        
                        }
                    //End base case
                
                }
                for(int m = 1; m<work_nb; m++)
                {
                    ix = index[m*2+0];
                    iy = index[m*2+1];
                    real a0 = work[BLOCKSIZE*BLOCKSIZE*m+0];
                    real a1 = work[BLOCKSIZE*BLOCKSIZE*m+1];
                    real a2 = work[BLOCKSIZE*BLOCKSIZE*m+2];
                    real a3 = work[BLOCKSIZE*BLOCKSIZE*m+3];
                    real a4 = work[BLOCKSIZE*BLOCKSIZE*m+4];
                    real a5 = work[BLOCKSIZE*BLOCKSIZE*m+5];
                    real a6 = work[BLOCKSIZE*BLOCKSIZE*m+6];
                    real a7 = work[BLOCKSIZE*BLOCKSIZE*m+7];
                    real a8 = work[BLOCKSIZE*BLOCKSIZE*m+8];
                    #if BLOCKSIZE ==4 
                    real a9 = work[BLOCKSIZE*BLOCKSIZE*m+9];
                    real aA = work[BLOCKSIZE*BLOCKSIZE*m+10];
                    real aB = work[BLOCKSIZE*BLOCKSIZE*m+11];
                    real aC = work[BLOCKSIZE*BLOCKSIZE*m+12];
                    real aD = work[BLOCKSIZE*BLOCKSIZE*m+13];
                    real aE = work[BLOCKSIZE*BLOCKSIZE*m+14];
                    real aF = work[BLOCKSIZE*BLOCKSIZE*m+15];
                    #endif
                    //std::cout<<"thread : "<<thNb<<" continuing "<<"( "<<ix<<","<<iy<<" )\n";
                    int delt_x = ix - ix_p;
                    int delt_y = iy - iy_p;
                    //std::cout<<"deltx: "<<delt_x<<" delt_y: "<<delt_y<<" \n";
                    //std::cout<<"thread "<<thNb<<" acumulators : "<<acc_00<<" "<<acc_01<<" "<<acc_02<<" "<<acc_03<<" "<<acc_10<<" "<<acc_11<<" "<<acc_12<<" "<<acc_13<<"\n";
                    if(delt_x !=0)
                        {
                            X2 = X + ix;
                            Y1_[0] += acc_00;
                            Y1_[1] += acc_01;
                            Y1_[2] += acc_02;
                            #if BLOCKSIZE == 4
                            Y1_[3] += acc_03;
                            acc_03 = 0;
                            #endif 
                            acc_00 = 0;
                            acc_01 = 0;
                            acc_02 = 0;
                            
                            Y1_= Y_base+ ix;
                            //std::cout<<"Y1_ set to "<<ix<<"\n";
                        }
                    if(delt_y !=0)
                        {
                           
                            X1 = X + iy;
                            Y2_[0] += acc_10;
                            Y2_[1] += acc_11;
                            Y2_[2] += acc_12;
                            #if BLOCKSIZE ==4
                            Y2_[3] += acc_13;
                            acc_13 = 0;
                            #endif 
                            acc_10 = 0;
                            acc_11 = 0;
                            acc_12 = 0;
                            
                            Y2_= Y_swap+iy;
                            //std::cout<<"Y2_ set to "<<iy<<"\n";
                        }
                        

                    if(ix-iy !=0)
                        {
                            
                            #if BLOCKSIZE==4
                            real y10 = 0;
                            real y11 = 0;
                            real y12=0;
                            real y13=0;
                            real y20 = 0;
                            real y21 = 0;
                            real y22 = 0;
                            real y23=0;
                            real r10 = X1[0];
                            real r11 = X1[1];
                            real r12 = X1[2];
                            real r13 = X1[3];
                            real r20 = X2[0];
                            real r21 = X2[1];
                            real r22 = X2[2];
                            real r23 = X2[3];
                            {real & c = a0;
                                
                                y10 += c*r10;
                                y20 += c*r20;
                            }
                            {
                                real & c = a1;
                                y11 += c * r10;
                                y20 += c * r21;
                            }
                            {
                                real & c = a4;
                                y10 += c * r11;
                                y21 += c * r20;
                            }
                            {
                                real & c = a2;
                                y12 += c * r10;
                                y20 += c * r22;
                            }
                            {
                                real& c = a8;
                                y10 += c * r12;
                                y22 += c* r20;
                            }
                            {
                                real & c = a3;
                                y13 += c * r10;
                                y20 += c * r23;
                            }
                            {
                                real & c  = aC;
                                y10 += c *r13;
                                y23 += c * r20;
                            }
                            {
                                real & c = a5;
                                y11+= c*r11;
                                y21 += c *r21;
                            }
                            {
                                real & c=  a6;
                                y12 += c*r11;
                                y21 += c*r22;
                            }
                            {
                                real & c = a9;
                                y11 += c *r12;
                                y22 += c *r21;
                            }
                            {
                                real &c  = a7;
                                y13 += c *r11;
                                y21 += c * r23;
                            }
                            {
                                real & c  = aD;
                                y11 +=c*r13;
                                y23 += c *r21;
                            }
                            {
                                real & c  = aA;
                                y12 += c * r12;
                                y22 += c* r22;
                            }
                            {
                                real & c = aE;
                                y12 += c * r13;
                                y23 += c * r22;
                            }
                            {
                                real & c = aB;
                                y13 += c * r12;
                                y22 += c *r23;
                            }
                            {
                                real & c = aF;
                                y13 += c * r13;
                                y23 += c *r23;
                            }
                           
                          
                            
                            acc_00 += y10;
                            acc_01 += y11;
                            acc_02 += y12;
                            acc_03 += y13;
                            acc_10 += y20;
                            acc_11 += y21;
                            acc_12 += y22;
                            acc_13 += y23;

                            #endif 
                            #if BLOCKSIZE==3
                                real y10 = 0;
                                real y11 = 0;
                                real y12=  0;
                                real y20 = 0;
                                real y21 = 0;
                                real y22 = 0;
                                real r10 = X1[0];
                                real r11 = X1[1];
                                real r12 = X1[2];
                                real r20 = X2[0];
                                real r21 = X2[1];
                                real r22 = X2[2];

                                    {real & c = a0;
                                        
                                        y10 += c*r10;
                                        y20 += c*r20;
                                    }
                                    {
                                        real & c = a1;
                                        y11 += c * r10;
                                        y20 += c * r21;
                                    }
                                    {
                                        real & c = a3;
                                        y10 += c * r11;
                                        y21 += c * r20;
                                    }
                                    {
                                        real & c = a2;
                                        y12 += c * r10;
                                        y20 += c * r22;
                                    }
                                    {
                                        real& c = a6;
                                        y10 += c * r12;
                                        y22 += c* r20;
                                    }
                                    {
                                        real & c = a4;
                                        y11 += c * r11;
                                        y21 += c * r21;
                                    }
                                    {
                                        real & c = a7;
                                        y11+= c*r12;
                                        y22 += c *r21;
                                    }
                                    {
                                        real & c = a5;
                                        y12 += c *r11;
                                        y21 += c *r22;
                                    }
                                    {
                                        real &c  = a8;
                                        y12 += c *r12;
                                        y22 += c * r22;
                                    }
                                acc_00 += y10;
                                acc_01 += y11;
                                acc_02 += y12;
                                acc_10 += y20;
                                acc_11 += y21;
                                acc_12 += y22;
                                  

                            #endif 
                            #if BLOCKSIZE==2
                            #endif 
                            

                       
                        }
                    
                    else
                        {
                            #if BLOCKSIZE ==4
                            acc_10 += a0 * X1[0] + a4 * X1[1] + a8 * X1[2] + aC * X1[3];
                            acc_11 += a1 * X1[0] + a5 * X1[1] + a9 * X1[2] + aD * X1[3];
                            acc_12 += a2 * X1[0] + a6 * X1[1] + aA * X1[2] + aE * X1[3];
                            acc_13 += a3 * X1[0] + a7 * X1[1] + aB * X1[2] + aF * X1[3];
                            #endif
                            #if BLOCKSIZE == 3
                            acc_10 += a0 * X1[0] + a3 * X1[1] + a6 * X1[2];
                            acc_11 += a1 * X1[0] + a4 * X1[1] + a7 * X1[2];
                            acc_12 += a2 * X1[0] + a5 * X1[1] + a8 * X1[2];
                            #endif 
                            
                        }
                ix_p = ix;
                iy_p = iy;
                }
                
               
                if(work_nb >=1)
                {
                Y1_[0] += acc_00;
                Y1_[1] += acc_01;
                Y1_[2] += acc_02;
                
                Y2_[0] += acc_10;
                Y2_[1] += acc_11;
                Y2_[2] += acc_12;
                #if BLOCKSIZE ==4
                    Y2_[3] += acc_13;
                    Y1_[3] += acc_03;
                #endif 
                }
                barrier2.arrive_and_wait();
                j++;
            }
        }

    void vecMulAddnTimes2(const real*X_calc, real*Y, int n_time)
        {

        X= X_calc;
        using milli = std::chrono::milliseconds;
        std::cout<<"[INFO]Starting MT calculation\n";
        auto start = std::chrono::high_resolution_clock::now();
        generateBlocks();
        auto stop = std::chrono::high_resolution_clock::now();
       

        std::cout<<"[INFO] Time spent in generateBlocks"<<std::chrono::duration_cast<milli>(stop - start).count()<<" ms\n";
        try
        {
            std::vector<std::thread> threads;
            std::barrier bar(nbThreads);
            std::barrier bar2(nbThreads);
            for(int i=0;i<nbThreads; i++)
            {
                threads.emplace_back(&MatSymBMtInstance::work2, this, std::ref(bar), std::ref(bar2), n_time);
            }
            
            
            
           
            
            for (auto& thread : threads)
            {
                thread.join();
            }
         
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception occurred: " << e.what() << std::endl;
        
        }
        for(int i=0;i<big_mat_size;i++)
        {
                
                Y[i]+= Y1[i]+Y2[i];
                
        }
        
    
        }

    void vecMulAddnTimes4(const real*X_calc, real*Y, int n_time)
        {

        X= X_calc;
        using milli = std::chrono::milliseconds;
        std::cout<<"[INFO]Starting MT calculation\n";
        auto start = std::chrono::high_resolution_clock::now();
        generateBlocks();
        auto stop = std::chrono::high_resolution_clock::now();
       

        std::cout<<"[INFO] Time spent in generateBlocks"<<std::chrono::duration_cast<milli>(stop - start).count()<<" ms\n";
        try
        {
            std::vector<std::thread> threads;
            std::barrier bar(nbThreads);
            std::barrier bar2(nbThreads);
            for(int i=0;i<nbThreads; i++)
            {
                threads.emplace_back(&MatSymBMtInstance::work4, this, std::ref(bar), std::ref(bar2), n_time);
            }
            
            
            
           
            
            for (auto& thread : threads)
            {
                thread.join();
            }
         
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception occurred: " << e.what() << std::endl;
        
        }
        for(int i=0;i<big_mat_size;i++)
        {
                
                Y[i]+= Y1[i]+Y2[i];
                
        }
        
    
        }




};

#endif /* MatSymBMtInstance_hpp */

