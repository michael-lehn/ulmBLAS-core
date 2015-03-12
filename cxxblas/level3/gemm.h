#ifndef CXXBLAS_LEVEL3_GEMM_H
#define CXXBLAS_LEVEL3_GEMM_H 1

#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB,
          typename Beta, typename TC>
    void
    gemm(bool         transA,
         bool         transB,
         IndexType    m,
         IndexType    n,
         IndexType    k,
         const Alpha  &alpha,
         bool         conjA,
         const TA     *A,
         IndexType    incRowA,
         IndexType    incColA,
         bool         conjB,
         const TB     *B,
         IndexType    incRowB,
         IndexType    incColB,
         const Beta   &beta,
         TC           *C,
         IndexType    incRowC,
         IndexType    incColC);

} // namespace cxxblas

#endif // CXXBLAS_LEVEL3_GEMM_H

#include <cxxblas/level3/gemm.tcc>
