//
//  MatSymBMtInstance.hpp
//  MatrixCalculation
//
//  Created by Moreau Enzo on 17/04/2024.
//

#ifndef MatSymBMtInstance_hpp
#define MatSymBMtInstance_hpp
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
    int blocksize = 4;
    real* Y2;
    
    std::ofstream tfile;//Debugging purposes only
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
    void generateBlocks2()
    {
        //tfile.open("test.txt", std::ios::app);
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
                //std::cout << phase_tab[i] << " ";
            }
            for(int i=1 ; i<nbThreads/2;i++)
            {
                phase_tab[nbThreads/2 + i] = nbThreads/2 - i;
                //std::cout << phase_tab[nbThreads/2 + i] << " ";
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
        threadBlockSize = (int) matrix->size() / (S_BLOCK_SIZE * nbThreads);
        int* phase_counter = (int*) malloc(sizeof(int)*phase_number);
        work_lengths = (int**) malloc(sizeof(int*)*phase_number);
        for(int i=0; i<phase_number;i++)
        {
            work_lengths[i] = (int*) malloc(sizeof(int)*nbThreads);
        }
        //Init
        
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
        
        int i =0;
        while(true)//Fetching blocks, assigned them to a phase each
        {
            int indX, indY;
            SparMatSymBlk::Column* col = &matrix->column_[i];
            if(col->nbb_>0)
            {
                //std::cout<<"Column number "<<matrix->colidx_[i]<<"->ind_YB : "<<(int) matrix->colidx_[i]/threadBlockSize<<"\n";
                indY = matrix->colidx_[i]/threadBlockSize;
                int indice_col = matrix->colidx_[i];
                
                for(int j=0; j<col->nbb_; j++)
                {
                    int indice_ligne = col->inx_[j];
                    //std::cout<<"Block position :"<<col->inx_[j]<<"-> ind_XB : "<<(int) col->inx_[j]/threadBlockSize<<"ind_YB:"<<indY<<"\n";
                    indX =col->inx_[j]/threadBlockSize;
                
                

                
                    int phase =phase_tab[indX-indY];
                   // std::cout<<"Chosen phase: "<<phase<<"\n";
                    //phasesTemp[phase][phase_counter[phase]%nbThreads]->push_front(Tuple2( indice_ligne, indice_col,swap, &matrix->block(col->inx_[j],matrix->colidx_[i])));
                    Block data;
                    Matrix44 matrix= col->blk_[j];
                    data.a0 = matrix.val[0x0];
                    data.a1 = matrix.val[0x1];
                    data.a2 = matrix.val[0x2];
                    data.a3 = matrix.val[0x3];
                    data.a4 = matrix.val[0x4];
                    data.a5 = matrix.val[0x5];
                    data.a6 = matrix.val[0x6];
                    data.a7 = matrix.val[0x7];
                    data.a8 = matrix.val[0x8];
                    data.a9 = matrix.val[0x9];
                    data.aA = matrix.val[0xA];
                    data.aB = matrix.val[0xB];
                    data.aC = matrix.val[0xC];
                    data.aD = matrix.val[0xD];
                    data.aE = matrix.val[0xE];
                    data.aF = matrix.val[0xF];
                    
                    data.index_x = indice_ligne * blocksize;
                    data.index_y = indice_col * blocksize;
                    
                  
                    phasesTemp[phase][indX]->push_front(data);
                  
                    
                    phase_counter[phase]++;
                    
                }
            }
            i = matrix->colidx_[i+1];
            if(i>=matrix->rsize_)
            {
                
                break;
            }
            
        }
        //std::ofstream outfile1;
        //outfile1.open("error.txt", std::ios::app);
        for(int i=0; i<phase_number; i++)//Putting in order for threads
        {
            for(int j=0; j<nbThreads;j++)
            {
                
                
                int t_work_length = (int)phasesTemp[i][j]->size();
                
                work_lengths[i][j] = t_work_length;
                
                workingPhases[i][j] = (real*) malloc(sizeof(real)*t_work_length*16);
                workingIndexes[i][j] = (int*) malloc(sizeof(int)*t_work_length*2);
                
                
                for(int m=0; m<t_work_length; m++)
                {
                    
                   
                    Block first= phasesTemp[i][j]->front();
                    phasesTemp[i][j]->pop_front();
                    
                    workingPhases[i][j][16*m] = first.a0;
                    workingPhases[i][j][16*m+1] = first.a1;
                    workingPhases[i][j][16*m+2] = first.a2;
                    workingPhases[i][j][16*m+3] = first.a3;
                    workingPhases[i][j][16*m+4] = first.a4;
                    workingPhases[i][j][16*m+5] = first.a5;
                    workingPhases[i][j][16*m+6] = first.a6;
                    workingPhases[i][j][16*m+7] = first.a7;
                    workingPhases[i][j][16*m+8] = first.a8;
                    workingPhases[i][j][16*m+9] = first.a9;
                    workingPhases[i][j][16*m+10] = first.aA;
                    workingPhases[i][j][16*m+11] = first.aB;
                    workingPhases[i][j][16*m+12] = first.aC;
                    workingPhases[i][j][16*m+13] = first.aD;
                    workingPhases[i][j][16*m+14] = first.aE;
                    workingPhases[i][j][16*m+15] = first.aF;
                    workingIndexes[i][j][2*m] = first.index_x;
                    workingIndexes[i][j][2*m+1] = first.index_y;
                    //outfile1<<"Block sent : ("<<first.index_x<<", "<<first.index_y<<") ("<<i<<","<<j<<")"<<"\n";
                    
                    
                }
                
                
            }
        }

        //outfile1<<"\n";
        //Must add a way to remove empty phases (happens a lot)
        for(int i =0; i<phase_number; i++)
        {
            //is the phase i Empty?
            int count = 0;
            
            for(int thNb = 0; thNb<nbThreads; thNb++)
            {
                count += work_lengths[i][thNb];//Count total number of work in the phase
            }
            if(count>0)
            {
                valid_phase_index.emplace_back(i);
            }
        }
        //We now know what phase we just do and which one we should not do

    }
    void generateBlocks3()
    {
        //tfile.open("test.txt", std::ios::app);
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
                //std::cout << phase_tab[i] << " ";
            }
            for(int i=1 ; i<nbThreads/2;i++)
            {
                phase_tab[nbThreads/2 + i] = nbThreads/2 - i;
                //std::cout << phase_tab[nbThreads/2 + i] << " ";
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
        threadBlockSize = (int) matrix->size() / (S_BLOCK_SIZE * nbThreads);
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
        
        int i =0;
        while(true)//Fetching blocks, assigned them to a phase each
        {
            int indX, indY;
            SparMatSymBlk::Column* col = &matrix->column_[i];
            if(col->nbb_>0)
            {
                //std::cout<<"Column number "<<matrix->colidx_[i]<<"->ind_YB : "<<(int) matrix->colidx_[i]/threadBlockSize<<"\n";
                indY = matrix->colidx_[i]/threadBlockSize;
                int indice_col = matrix->colidx_[i];
                
                for(int j=0; j<col->nbb_; j++)
                {
                    int indice_ligne = col->inx_[j];
                    //std::cout<<"Block position :"<<col->inx_[j]<<"-> ind_XB : "<<(int) col->inx_[j]/threadBlockSize<<"ind_YB:"<<indY<<"\n";
                    indX =col->inx_[j]/threadBlockSize;
                
                

                
                    int phase =phase_tab[indX-indY];
                   // std::cout<<"Chosen phase: "<<phase<<"\n";
                    //phasesTemp[phase][phase_counter[phase]%nbThreads]->push_front(Tuple2( indice_ligne, indice_col,swap, &matrix->block(col->inx_[j],matrix->colidx_[i])));
                    Block data;
                    Matrix44 matrix= col->blk_[j];
                    data.a0 = matrix.val[0x0];
                    data.a1 = matrix.val[0x1];
                    data.a2 = matrix.val[0x2];
                    data.a3 = matrix.val[0x3];
                    data.a4 = matrix.val[0x4];
                    data.a5 = matrix.val[0x5];
                    data.a6 = matrix.val[0x6];
                    data.a7 = matrix.val[0x7];
                    data.a8 = matrix.val[0x8];
                    data.a9 = matrix.val[0x9];
                    data.aA = matrix.val[0xA];
                    data.aB = matrix.val[0xB];
                    data.aC = matrix.val[0xC];
                    data.aD = matrix.val[0xD];
                    data.aE = matrix.val[0xE];
                    data.aF = matrix.val[0xF];
                    
                    data.index_x = indice_ligne * blocksize;
                    data.index_y = indice_col * blocksize;
                    
                  
                    phasesTemp[phase][indX]->push_front(data);
                  
                    
                    phase_counter[phase]++;
                    
                }
            }
            i = matrix->colidx_[i+1];
            if(i>=matrix->rsize_)
            {
                
                break;
            }
            
        }
        
        newPhases = (real**) malloc(phase_number*sizeof(real*));
        newIndexes = (int**) malloc(phase_number*sizeof(int*));
        for(int i=0; i<phase_number; i++)//Putting in order for threads
        {
            int maxCount =0;
            for(int j=0; j<nbThreads;j++)
            {
                
                
                int t_work_length = (int)phasesTemp[i][j]->size();
                if(t_work_length > maxCount)
                {
                    maxCount = t_work_length;
                }
                work_lengths[i][j] = t_work_length;
            }
            if(maxCount ==0)
            {
                newPhases[i] = (real*) malloc(1*sizeof(real));
                newIndexes[i] = (int*) malloc(1*sizeof(int));
            }
            else
            {
                //std::cout<<"Phase:"<<i<<"Max is "<<maxCount*2*nbThreads<<"\n";
                newPhases[i] = (real*) malloc(maxCount*16*nbThreads*sizeof(real));
                newIndexes[i] = (int*) malloc(maxCount*2*nbThreads*sizeof(int));
                for(int j=0; j<nbThreads;j++)
                {
                    for(int m=0; m<work_lengths[i][j];m++)
                    {
                    Block block = phasesTemp[i][j]->front();
                    //std::cout<<"Send block at "
                    phasesTemp[i][j]->pop_front();

                    newPhases[i][16*(m*nbThreads+j)] = block.a0;
                    newPhases[i][16*(m*nbThreads+j)+1] = block.a1;
                    newPhases[i][16*(m*nbThreads+j)+2] = block.a2;
                    newPhases[i][16*(m*nbThreads+j)+3] = block.a3;
                    newPhases[i][16*(m*nbThreads+j)+4] = block.a4;
                    newPhases[i][16*(m*nbThreads+j)+5] = block.a5;
                    newPhases[i][16*(m*nbThreads+j)+6] = block.a6;
                    newPhases[i][16*(m*nbThreads+j)+7] = block.a7;
                    newPhases[i][16*(m*nbThreads+j)+8] = block.a8;
                    newPhases[i][16*(m*nbThreads+j)+9] = block.a9;
                    newPhases[i][16*(m*nbThreads+j)+10] = block.aA;
                    newPhases[i][16*(m*nbThreads+j)+11] = block.aB;
                    newPhases[i][16*(m*nbThreads+j)+12] = block.aC;
                    newPhases[i][16*(m*nbThreads+j)+13] = block.aD;
                    newPhases[i][16*(m*nbThreads+j)+14] = block.aE;
                    newPhases[i][16*(m*nbThreads+j)+15] = block.aF;
                    newIndexes[i][2*(m*nbThreads+j)+0]= block.index_x;
                    newIndexes[i][2*(m*nbThreads+j)+1]= block.index_y;
                    }
                }
            }
             
        }
        for(int i =0; i<phase_number; i++)
        {
            //is the phase i Empty?
            int count = 0;
            
            for(int thNb = 0; thNb<nbThreads; thNb++)
            {
                count += work_lengths[i][thNb];//Count total number of work in the phase
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
        //std::cout<<"Thread number:  "<<thread_nb<<"\n";
        _mutex.unlock();
        
        for(int m = 0; m<n_work;m++)
        {
            barrier.arrive_and_wait();
            workThread2(barrier2, thread_nb);
            barrier.arrive_and_wait();
        }
   
      
        
        
        
        
    } 
    void work3(std::barrier<>& barrier,std::barrier<>& barrier2, int n_work)
    {
       
        _mutex.lock();
        int thread_nb = thread_number();
   
        _mutex.unlock();
        
        for(int m = 0; m<n_work;m++)
        {
            
            workThread3(barrier2, thread_nb);
            barrier.arrive_and_wait();
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
        int time1 = 0;
        int time2 = 0;
        int nb1 = 0;
        int nb2 = 0;
        
        while(j<valid_phase_index.size())
        {
            
            

            int k  = valid_phase_index[j];
            work =  workingPhases[k][thNb];
            
            work_nb = work_lengths[k][thNb];
            int* index = workingIndexes[k][thNb];
            
            
            for(int m = 0; m<work_nb; m++)
            {
                auto start = std::chrono::high_resolution_clock::now();
                
                int ix = index[m*2+0];
                int iy = index[m*2+1];

        
               
              
                real a0 = work[16*m+0];
                real a1 = work[16*m+1];
                real a2 = work[16*m+2];
                real a3 = work[16*m+3];
                real a4 = work[16*m+4];
                real a5 = work[16*m+5];
                real a6 = work[16*m+6];
                real a7 = work[16*m+7];
                real a8 = work[16*m+8];
                real a9 = work[16*m+9];
                real aA = work[16*m+10];
                real aB = work[16*m+11];
                real aC = work[16*m+12];
                real aD = work[16*m+13];
                real aE = work[16*m+14];
                real aF = work[16*m+15];
                
                if(ix-iy <= threadBlockSize || ix-iy==0)
                {
                    X1 = X  + iy; //no swap
                    X2 = X  + ix;
                    Y1_= Y1 + ix;
                    Y2_= Y2 + iy;
                    
                }
                else
                {
                    //swap
                    X1 = X  + iy;
                    X2 = X  + ix;
                    Y1_= Y2 + ix;
                    Y2_= Y1 + iy;
                    
                    
                }
                
                if(ix-iy !=0)
                {
                    
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
                    
                 
                    

                  
                }
                else
                {
                    
                    Y2_[0] += a0 * X1[0] + a4 * X1[1] + a8 * X1[2] + aC * X1[3];
                    Y2_[1] += a1 * X1[0] + a5 * X1[1] + a9 * X1[2] + aD * X1[3];
                    Y2_[2] += a2 * X1[0] + a6 * X1[1] + aA * X1[2] + aE * X1[3];
                    Y2_[3] += a3 * X1[0] + a7 * X1[1] + aB * X1[2] + aF * X1[3];
                  
                }
            }
            
           
            j++;
            barrier2.arrive_and_wait();
            
        }
}
    void workThread3(std::barrier<>& barrier2, int thNb)
    {
        

        int active_work_count;
        int j = 0;
        
        real* work;
        int work_nb;
        
        const real* X1;
        const real* X2;
        real* Y1_;
        real* Y2_;
        int time1 = 0;
        int time2 = 0;
        int nb1 = 0;
        int nb2 = 0;
        //std::cout<<"Checkpoint1";
        while(j<valid_phase_index.size())
        {
            
            
         
           
            int k  = valid_phase_index[j];
            work =  newPhases[k];//Fetch the phase
            
            work_nb = work_lengths[k][thNb];
            int* index = newIndexes[k];
            
           
            
            for(int m = 0; m<work_nb; m++)
            {
            
                
                
                int ix = index[(thNb+m*nbThreads)*2+0];
                int iy = index[(thNb+m*nbThreads)*2+1];
              
                
                
                
                
        
               
               
                real a0 = work[16*(m*nbThreads+thNb)+0];
                real a1 = work[16*(m*nbThreads+thNb)+1];
                real a2 = work[16*(m*nbThreads+thNb)+2];
                real a3 = work[16*(m*nbThreads+thNb)+3];
                real a4 = work[16*(m*nbThreads+thNb)+4];
                real a5 = work[16*(m*nbThreads+thNb)+5];
                real a6 = work[16*(m*nbThreads+thNb)+6];
                real a7 = work[16*(m*nbThreads+thNb)+7];
                real a8 = work[16*(m*nbThreads+thNb)+8];
                real a9 = work[16*(m*nbThreads+thNb)+9];
                real aA = work[16*(m*nbThreads+thNb)+10];
                real aB = work[16*(m*nbThreads+thNb)+11];
                real aC = work[16*(m*nbThreads+thNb)+12];
                real aD = work[16*(m*nbThreads+thNb)+13];
                real aE = work[16*(m*nbThreads+thNb)+14];
                real aF = work[16*(m*nbThreads+thNb)+15];
                
                if(ix-iy <= threadBlockSize || ix-iy==0)
                {
                    X1 = X  + iy; //no swap
                    X2 = X  + ix;
                    Y1_= Y1 + ix;
                    Y2_= Y2 + iy;
                    
                }
                else
                {
                    //swap
                    X1 = X  + iy;
                    X2 = X  + ix;
                    Y1_= Y2 + ix;
                    Y2_= Y1 + iy;
                    
                    
                }
                
                if(ix-iy !=0)
                {
                    
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
                    
                 
                    

                   
                }
                else
                {
                    Y2_[0] += a0 * X1[0] + a4 * X1[1] + a8 * X1[2] + aC * X1[3];
                    Y2_[1] += a1 * X1[0] + a5 * X1[1] + a9 * X1[2] + aD * X1[3];
                    Y2_[2] += a2 * X1[0] + a6 * X1[1] + aA * X1[2] + aE * X1[3];
                    Y2_[3] += a3 * X1[0] + a7 * X1[1] + aB * X1[2] + aF * X1[3];
 
                
                }
            }
            
           
            j++;
            barrier2.arrive_and_wait();
            
            
        }
}
   
    void vecMulAddnTimes2(const real*X_calc, real*Y, int n_time)
    {

        X= X_calc;
        using milli = std::chrono::milliseconds;
        std::cout<<"[INFO]Starting MT calculation\n";
        auto start = std::chrono::high_resolution_clock::now();
        generateBlocks2();
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
            // Handle exception here if needed
            
        }
        for(int i=0;i<big_mat_size;i++)
        {
                
                Y[i]+= Y1[i]+Y2[i];
                
        }
        
    
    }
    void vecMulAddnTimes3(const real*X_calc, real*Y, int n_time)
    {

        X= X_calc;
        using milli = std::chrono::milliseconds;
        std::cout<<"[INFO]Starting MT calculation\n";
        auto start = std::chrono::high_resolution_clock::now();
        generateBlocks3();
        auto stop = std::chrono::high_resolution_clock::now();
       

        std::cout<<"[INFO] Time spent in generateBlocks "<<std::chrono::duration_cast<milli>(stop - start).count()<<" ms\n";
        try
        {
            std::vector<std::thread> threads;
            std::barrier bar(nbThreads);
            std::barrier bar2(nbThreads);
            for(int i=0;i<nbThreads; i++)
            {
                threads.emplace_back(&MatSymBMtInstance::work3, this, std::ref(bar), std::ref(bar2), n_time);
            }
            
            
            
           
            
            for (auto& thread : threads)
            {
                thread.join();
            }
         
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception occurred: " << e.what() << std::endl;
            // Handle exception here if needed
            
        }
        for(int i=0;i<big_mat_size;i++)
        {
                
                Y[i]+= Y1[i]+Y2[i];
                
        }
        
    
    }
    

};

#endif /* MatSymBMtInstance_hpp */


