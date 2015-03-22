#ifndef ULMBLAS_CXXBLAS_LEVEL3_HERK_H
#define ULMBLAS_CXXBLAS_LEVEL3_HERK_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename Beta,
          typename TC>
    void
    herk(IndexType    n,
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

#endif // ULMBLAS_CXXBLAS_LEVEL3_HERK_H

#include <ulmblas/cxxblas/level3/herk.tcc>
