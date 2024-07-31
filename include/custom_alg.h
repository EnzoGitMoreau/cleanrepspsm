#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#ifdef RSB
        #include <rsb.hpp>
    #endif

std::vector<std::pair<int, int>> select_random_points(int n, int k) {
    // Create a vector of pairs representing all points in the n x n grid
    std::vector<std::pair<int, int>> all_points;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <=i; j++) {
            all_points.push_back(std::make_pair(i, j));
        }
    }
    k = k/2;//We are in symmetric representation -> 20% of full block = 20/2 % of possible coordinates  


    // Shuffle the vector of points
    std::mt19937 rng(std::random_device{}());
    std::shuffle(all_points.begin(), all_points.end(), rng);

    
    std::vector<std::pair<int, int>> selected_points(all_points.begin(), all_points.begin() + k);
    return selected_points;
}

void add_block_to_pos(double* values, std::vector<std::pair<int, int>> pairs, int lda) {
    for (size_t k = 0; k < pairs.size(); k++) { // Iterate over all the pairs in the vector
        int row = pairs[k].first;
        int col = pairs[k].second;
        for (int i = 0; i < BLOCKSIZE; i++) { // Iterate over the 4x4 block
            for (int j = 0; j < BLOCKSIZE; j++) 
			{
                values[(row*BLOCKSIZE + i) * lda + (col*BLOCKSIZE + j)] = (i + j) / 4.0;
				values[(col*BLOCKSIZE + j) * lda + row*BLOCKSIZE+i] = (i + j) / 4.0;
            }
        }
    }
}

void add_block_to_pos_std(SparMatSymBlk* matrix, std::vector<std::pair<int, int>> pairs, int lda)
{
	for (size_t k = 0; k < pairs.size(); k++) 
	{ 
        int row = pairs[k].first;
        int col = pairs[k].second;
        #if BLOCKSIZE == 4
        Matrix44* block= new Matrix44(1,2/4.0,3/4.0,4/4.0,1/4.0,2/4.0,3/4.0,4/4,1/4.0,2/4.0,3/4.0,4/4.0,1/4.0,2/4.0,3/4.0,4/4.0);
        #endif 
        #if BLOCKSIZE == 3 
        Matrix33* block= new Matrix33(1,2/4.0,3/4.0,4/4.0,1/4.0,2/4.0,3/4.0,4/4,1/4.0);
        #endif
        #if BLOCKSIZE ==2 
        Matrix22* block= new Matrix22(1,2/4.0,3/4.0,4/4.0);
        #endif
		matrix->block(row,col).add_full(*block);
    }

}
#ifdef RSB
void add_block_to_pos_rsb(rsb::RsbMatrix<double>** mtx_ptr, SparMatSymBlk* matrix, std::vector<std::pair<int, int>> pairs, int lda)
{
    const rsb_coo_idx_t nrA { lda}, ncA { lda };
	*mtx_ptr = new rsb::RsbMatrix<double>(nrA, ncA);
	for (size_t k = 0; k < pairs.size(); k++) 
	{ 
        int row = pairs[k].first;
        int col = pairs[k].second;
        #if BLOCKSIZE == 4
        Matrix44* block= new Matrix44(1,2/4.0,3/4.0,4/4.0,1/4.0,2/4.0,3/4.0,4/4,1/4.0,2/4.0,3/4.0,4/4.0,1/4.0,2/4.0,3/4.0,4/4.0);
        #endif 
        #if BLOCKSIZE == 3 
        Matrix33* block= new Matrix33(1,2,3,4,5,6,7,8,9);//to change
        #endif
        #if BLOCKSIZE ==2 
        Matrix22* block= new Matrix22(1,2/4.0,3/4.0,4/4.0);
        #endif
		matrix->block(row,col).add_full(*block);
        for(size_t m = 0; m<BLOCKSIZE; m++)
        {
            for(size_t n=0; n<BLOCKSIZE; n++)
            {
                int i = row +m;
                int j = col + n;
            double errval = (*mtx_ptr)->get_val(i,j);
            if(!errval)
            {(*mtx_ptr)->set_val((m+n)/4.0, i,j);
            if(i!=j)
                {
                (*mtx_ptr)->set_val((m+n)/4.0, j,i);
                }}
            }
            
        }
        
        
    }
    (*mtx_ptr)->close();

}
#endif

#include <matsym.h>
void SPSMtoSym(MatrixSymmetric* dest, SparMatSymBlk* origin)
{
    
    for ( size_t j = origin->colidx_[0]; j <origin->rsize_; j = origin->colidx_[j+1] )
    {

        SparMatSymBlk::Column* col = &origin->column_[j];;
        for ( size_t n = 0; n < col->nbb_; ++n )
    {
        
        #if BLOCKSIZE ==3
        
        Matrix33 &M = col->blk_[n];
 
        int x = n*BLOCKSIZE;
        int y = j*BLOCKSIZE;
     

        dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[0]; 
        x++;
        dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[3]; 
        x++;
        dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[6]; 
        y++;

        x = n*BLOCKSIZE;
        dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[1]; 
        x++;
         dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[4]; 
        x++;
         dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[7]; 
        y++;
        
        x = n*BLOCKSIZE;
         dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[2]; 
        x++;
         dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[5]; 
        x++;
         dest->val[ std::max(x,y) + dest->size() * std::min(x,y) ] = M.val[8]; 
        y++;

       
        




        #endif 
        #if BLOCKSIZE ==4
        Matrix44 &M=col->blk_[n];
        #endif 
        #if BLOCKSIZE ==2
        Matrix22 & M = col.blk_[n];
        #endif
        //std::cout<< found block
    }

    }
}
