#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
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
        for (int i = 0; i < 4; i++) { // Iterate over the 4x4 block
            for (int j = 0; j < 4; j++) 
			{
                values[(row*4 + i) * lda + (col*4 + j)] = (i + j) / 4.0;
				values[(col*4 + j) * lda + row*4+i] = (i + j) / 4.0;
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
        Matrix44* block= new Matrix44(1,2/4.0,3/4.0,4/4.0,1/4.0,2/4.0,3/4.0,4/4,1/4.0,2/4.0,3/4.0,4/4.0,1/4.0,2/4.0,3/4.0,4/4.0);
		matrix->block(row,col).add_full(*block);
    }

}