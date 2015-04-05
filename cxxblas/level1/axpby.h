#ifndef ULMBLAS_CXXBLAS_LEVEL1_AXPBY_H
#define ULMBLAS_CXXBLAS_LEVEL1_AXPBY_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename Beta,
          typename TY>
    void
    axpby(IndexType    n,
          const Alpha  &alpha,
          const TX     *x,
          IndexType    incX,
          const Beta   &beta,
          TY           *y,
          IndexType    incY);

template <typename IndexType, typename Alpha, typename TX, typename Beta,
          typename TY>
    void
    acxpby(IndexType    n,
           const Alpha  &alpha,
           const TX     *x,
           IndexType    incX,
           const Beta   &beta,
           TY           *y,
           IndexType    incY);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/axpby.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_AXPBY_H
