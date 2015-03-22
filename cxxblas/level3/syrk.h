#ifndef ULMBLAS_CXXBLAS_LEVEL3_SYRK_H
#define ULMBLAS_CXXBLAS_LEVEL3_SYRK_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename Beta,
          typename TC>
    void
    syrk(IndexType    n,
         IndexType    k,
         const Alpha  &alpha,
         bool         transA,
         const TA     *A,
         IndexType    incRowA,
         IndexType    incColA,
         const Beta   &beta,
         bool         lowerC,
         TC           *C,
         IndexType    incRowC,
         IndexType    incColC);

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_SYRK_H

#include <ulmblas/cxxblas/level3/syrk.tcc>
