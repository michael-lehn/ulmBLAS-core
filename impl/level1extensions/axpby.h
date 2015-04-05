#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPBY_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPBY_H 1

namespace ulmBLAS {

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

} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/axpby.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPBY_H 1
