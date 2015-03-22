#ifndef ULMBLAS_CXXBLAS_LEVEL2_SPMV_H
#define ULMBLAS_CXXBLAS_LEVEL2_SPMV_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
    void
    spmv(IndexType    n,
         const Alpha  &alpha,
         bool         colMajorA,
         bool         lowerA,
         const TA     *AP,
         const TX     *x,
         IndexType    incX,
         const Beta   &beta,
         TY           *y,
         IndexType    incY);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/spmv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_SPMV_H
