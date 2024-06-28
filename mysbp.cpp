#include <rsb.hpp>
#include <vector>
#include <array>

int main() {
  rsb::RsbLib rsblib;
  const int nnzA { 7 }, nrhs { 2 };
  const int nrA { 6 }, ncA { 6 };
  const std::vector<int>    IA {0,1,2,3,4,5,1};
  const int                 JA [] = {0,1,2,3,4,5,0};
  const std::vector<double> VA {1,1,1,1,1,1,2}, B(nrhs * ncA,1);
  std::array<double,nrhs * nrA> C;
  const double alpha {2}, beta {1};

  // The three arrays IA, JA, VA form a COO (Coordinate) representation of a 6x6 matrix
  rsb::RsbMatrix<double> mtx(IA,JA,VA,nnzA); // Declarations of IA,JA,VA are all accepted via <span>

  mtx.spmm(RSB_TRANSPOSITION_N, alpha, nrhs, RSB_FLAG_WANT_ROW_MAJOR_ORDER, B, beta, C);
}