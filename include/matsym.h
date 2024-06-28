// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#ifndef MATSYM_H
#define MATSYM_H
#define NB_THREADS 2

//Nb thread = 6 doesnt work
//Nb thread = 4 doesnt work
#include "real.h"

#include <cstdio>
#ifdef MACOS 
#include <arm_neon.h>
#endif
#include <string>


class MatrixSymmetric final
{
private:
    
    /// leading dimension of array
    size_t dimension_;
    
    /// size of matrix
    size_t size_;

    /// size of memory which has been allocated
    size_t allocated_;
    
    // full upper triangle:
    real* val;
    
    // if 'false', destructor will not call delete[] val;
    bool in_charge;
    
public:
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSymmetric();
    
    
    /// constructor from an existing array
    MatrixSymmetric(size_t s)
    {
        val = nullptr;
        resize(s);
        dimension_ = s;
        val = new_real(s*s);
        zero_real(s*s, val);
        in_charge = true;
    }

    /// constructor from an existing array
    MatrixSymmetric(size_t s, real* array, size_t ldd)
    {
        free_real(val);
        size_ = s;
        dimension_ = ldd;
        val = array;
        in_charge = false;
    }
    
    /// default destructor

    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t alc);
    
    /// returns address of data array
    real* data() const { return val; }

    /// returns the address of element at (x, y), no allocation is done
    real* addr(size_t x, size_t y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(size_t i, size_t j);
    
    /// scale the matrix by a scalar factor
    void scale(real a);
    
    /// multiplication of a vector: Y = Y + M * X, dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y) const;
    
    /// 2D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso2D(const real* X, real* Y) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso3D(const real* X, real* Y) const;
    void vecMulAddBlock(const real* X, real*Y, int index_x, int index_j, int blocksize, int matsize) const;
    void transVecMulAddBlock(const real* X, real*Y, int index_x, int index_j, int blocksize, int matsize) const;
    void vecMulAddBlock2(const real * __restrict__ X, real *__restrict__ Y, int index_x, int index_j, int blocksize, int matsize) const;
    void transVecMulAddBlock2(const real* X, real*Y1, real*Y2, int index_x1, int index_y1, int blocksize, int matsize) const;
    void transVecMulAddBlock3(const real* __restrict__ X1, const real* __restrict__ X2,real* __restrict__ Y1,  real* __restrict__ Y2, int blocksize, int matsize, const real* __restrict__ valptr )const;
    /// true if matrix is non-zero
    void transVecMulAddBlock4(const real* __restrict__ X1, const real* __restrict__ X2,real* __restrict__ Y1,  real* __restrict__ Y2, int blocksize, int matsize, const real* __restrict__ valptr )const;
   #ifdef MACOS
    void matrix_multiply_4x4_neon(const real* __restrict__ X1, const real* __restrict__ X2,real* __restrict__ Y1,  real* __restrict__ Y2, int blocksize, int matsize, const real* __restrict__ valptr)const;
    void matrix_multiply_4x4_neonMid(const real* __restrict__ X1,real* __restrict__ Y1 ,int blocksize, int matsize, const real* __restrict__ valptr)const;
    
    void matrix_multiply_4x4_neon2(const real* __restrict__ X1, const real* __restrict__ X2,real* __restrict__ Y1,  real* __restrict__ Y2, int blocksize, int matsize, const real* __restrict__ valptr) const{
                // these are the columns A
            
                //Matrice 4x4
                
                float64x2_t A01;
                float64x2_t A02;
                float64x2_t A11;
                float64x2_t A12;
                float64x2_t A21;
                float64x2_t A22;
                float64x2_t A31;
                float64x2_t A32;
                float64x2_t X01;
                float64x2_t X02;
                float64x2_t X11;
                float64x2_t X12;
                float64x2_t Y21;
                float64x2_t Y11;
                float64x2_t Y12;
                float64x2_t Y22;
                
        float64x2_t partialSum;
        float64x2_t partialSum2;
        float64x2_t accumulator1;
        float64x2_t accumulator2;
                X01 = vld1q_f64(X1);
                X02 = vld1q_f64(X1+2);
                X11 = vld1q_f64(X2);
                X12 = vld1q_f64(X2+2);
       

                
                Y11 = vld1q_f64(Y1);
                Y12  = vld1q_f64(Y1+2);
       
        //accumulator1 = {0,0};
        //accumulator2 = {0,0};
                accumulator1 =vld1q_f64(Y1);
                accumulator2 = vld1q_f64(Y1+2);
        
                Y21 = vld1q_f64(Y2);
                Y22 = vld1q_f64(Y2+2);
               
        
                A01 = vld1q_f64(valptr);
                accumulator1 =vfmaq_n_f64(accumulator1, A01, vgetq_lane_f64(X01,0));
                partialSum = vmulq_f64(A01, X11);
                
                
                
                A02 = vld1q_f64(valptr+2);
                
                accumulator2 =vfmaq_n_f64(accumulator2, A02, vgetq_lane_f64(X01,0));
                partialSum = vfmaq_f64(partialSum, A02, X12);
        
        
        
                A11 = vld1q_f64(valptr+matsize);
                accumulator1 =vfmaq_n_f64(accumulator1, A11, vgetq_lane_f64(X01,1));
                partialSum2 =vmulq_f64(A11, X11);
        
        
        
                A12 = vld1q_f64(valptr+matsize+2);
                accumulator2 =vfmaq_n_f64(accumulator2, A12, vgetq_lane_f64(X01,1));
                partialSum2 = vfmaq_f64(partialSum2, A12, X12);
                ///Finishing calcul
                Y21 =vzip1q_f64(partialSum, partialSum2);
                partialSum = vzip2q_f64(partialSum, partialSum2);
                Y21 = vaddq_f64(partialSum, Y21);
                Y21 = vaddq_f64(vld1q_f64( Y2), Y21);
                vst1q_f64(Y2, Y21);
        
                A21 = vld1q_f64(valptr+2*matsize);
                accumulator1 =vfmaq_n_f64(accumulator1, A21, vgetq_lane_f64(X02,0));
                partialSum = vmulq_f64( A21, X11);
        
                A22 = vld1q_f64(valptr+2*matsize+2);
                accumulator2 =vfmaq_n_f64(accumulator2, A22, vgetq_lane_f64(X02,0));
                partialSum = vfmaq_f64(partialSum, A22, X12);
                
        
        
                A31 = vld1q_f64(valptr+3*matsize);
                accumulator1 =vfmaq_n_f64(accumulator1, A31, vgetq_lane_f64(X02,1));
                partialSum2 =vmulq_f64(A31, X11);
        
                A32 = vld1q_f64(valptr+3*matsize+2);
                accumulator2 =vfmaq_n_f64(accumulator2, A32, vgetq_lane_f64(X02,1));
                partialSum2 = vfmaq_f64(partialSum2, A32, X12);
                //Finishing calcul
                Y22 =vzip1q_f64(partialSum, partialSum2);
                partialSum = vzip2q_f64(partialSum, partialSum2);
                Y22 = vaddq_f64(partialSum, Y22);
                Y22 = vaddq_f64(vld1q_f64( Y2+2), Y22);
        vst1q_f64(Y2+2, Y22);
        vst1q_f64(Y1,accumulator1);
        vst1q_f64(Y1+2,accumulator2);
                
                
            
        }
   #endif
    void vecMulPerBlock( const real* X, real* Y,size_t blocksize, int nbThread ) const;
    
    
    bool notZero() const;
    
    /// number of element which are non-zero
    size_t nbElements(size_t start, size_t stop) const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
};

#endif
