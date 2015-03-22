#ifndef ULMBLAS_CXXBLAS_LEVEL2_SBMV_H
#define ULMBLAS_CXXBLAS_LEVEL2_SBMV_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
    void
    sbmv(IndexType    n,
         IndexType    k,
         const Alpha  &alpha,
         bool         colMajorA,
         bool         lowerA,
         const TA     *A,
         IndexType    ldA,
         const TX     *x,
         IndexType    incX,
         const Beta   &beta,
         TY           *y,
         IndexType    incY);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/sbmv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_SBMV_H

