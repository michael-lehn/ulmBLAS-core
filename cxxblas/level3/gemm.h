#ifndef ULMBLAS_CXXBLAS_LEVEL3_GEMM_H
#define ULMBLAS_CXXBLAS_LEVEL3_GEMM_H 1

#ifdef USE_EXTERNAL_BLAS
#include <ulmblas/external/blis/gemm.h>
#endif // USE_EXTERNAL_BLAS

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
         IndexType    incColC);

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_GEMM_H

#include <ulmblas/cxxblas/level3/gemm.tcc>
