#ifndef ULMBLAS_CXXBLAS_LEVEL3_SYR2K_H
#define ULMBLAS_CXXBLAS_LEVEL3_SYR2K_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB,
          typename Beta, typename TC>
    void
    syr2k(IndexType    n,
          IndexType    k,
          const Alpha  &alpha,
          bool         transAB,
          const TA     *A,
          IndexType    incRowA,
          IndexType    incColA,
          const TB     *B,
          IndexType    incRowB,
          IndexType    incColB,
          const Beta   &beta,
          bool         lowerC,
          TC           *C,
          IndexType    incRowC,
          IndexType    incColC);

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_SYR2K_H

#include <ulmblas/cxxblas/level3/syr2k.tcc>
