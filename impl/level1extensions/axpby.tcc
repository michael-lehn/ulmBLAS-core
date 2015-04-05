#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPBY_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPBY_TCC 1

#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/level1extensions/axpby.h>

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
      IndexType    incY)
{
    const IndexType    UnitStride(1);

    if (n<=0) {
        return;
    }

    if (incX==UnitStride && incY==UnitStride) {
        for (IndexType i=0; i<n; ++i) {
            y[i] = beta*y[i] + alpha*x[i];
        }
    } else {
        for (IndexType i=0; i<n; ++i) {
            y[i*incY] = beta*y[i*incY] + alpha*x[i*incX];
        }
    }
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
    const IndexType    UnitStride(1);

    if (n<=0) {
        return;
    }

    if (incX==UnitStride && incY==UnitStride) {
        for (IndexType i=0; i<n; ++i) {
            y[i] = beta*y[i] + alpha*x[i];
        }
    } else {
        for (IndexType i=0; i<n; ++i) {
            y[i*incY] = beta*y[i*incY] + alpha*conjugate(x[i*incX]);
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPBY_TCC 1
