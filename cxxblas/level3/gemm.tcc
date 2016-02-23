#ifndef ULMBLAS_CXXBLAS_LEVEL3_GEMM_TCC
#define ULMBLAS_CXXBLAS_LEVEL3_GEMM_TCC 1

#include <ulmblas/cxxblas/level3/gemm.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB,
          typename Beta, typename TC>
void
gemm(IndexType    m,
     IndexType    n,
     IndexType    k,
     const Alpha  &alpha,
     bool         transA,
     bool         conjA,
     const TA     *A,
     IndexType    incRowA_,
     IndexType    incColA_,
     bool         transB,
     bool         conjB,
     const TB     *B,
     IndexType    incRowB_,
     IndexType    incColB_,
     const Beta   &beta,
     TC           *C,
     IndexType    incRowC,
     IndexType    incColC)
{
    const IndexType incRowA = (!transA) ? incRowA_ : incColA_;
    const IndexType incColA = (!transA) ? incColA_ : incRowA_;

    const IndexType incRowB = (!transB) ? incRowB_ : incColB_;
    const IndexType incColB = (!transB) ? incColB_ : incRowB_;

    ulmBLAS::gemm(m, n, k,
                  alpha,
                  conjA, A, incRowA, incColA,
                  conjB, B, incRowB, incColB,
                  beta,
                  C, incRowC, incColC);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_GEMM_TCC
