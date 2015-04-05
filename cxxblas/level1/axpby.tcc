#ifndef ULMBLAS_CXXBLAS_LEVEL1_AXPBY_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_AXPBY_TCC 1

#include <ulmblas/cxxblas/level1/axpby.tcc>
#include <ulmblas/impl/level1extensions/axpby.h>

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
      IndexType    incY)
{
    ulmBLAS::axpby(n, alpha, x, incX, beta, y, incY);
}

template <typename IndexType, typename Alpha, typename TX, typename Beta,
          typename TY>
void
acxpby(IndexType    n,
       const Alpha  &alpha,
       const TX     *x,
       IndexType    incX,
       const Beta   &beta,
       TY           *y,
       IndexType    incY)
{
    ulmBLAS::acxpby(n, alpha, x, incX, beta, y, incY);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_AXPBY_TCC
