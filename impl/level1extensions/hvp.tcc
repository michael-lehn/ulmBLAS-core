#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_HVP_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_HVP_TCC 1

#include <ulmblas/impl/level1extensions/hvp.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename VX, typename VY>
void
ihvp(IndexType      n,
     const Alpha    &alpha,
     const VX       *x,
     IndexType      incX,
     VY             *y,
     IndexType      incY)
{
    const IndexType  UnitStride(1);

    if (incX==UnitStride && incY==UnitStride) {
        for (IndexType i=0; i<n; ++i) {
            y[i] *= alpha*x[i];
        }
    } else {
        for (IndexType i=0; i<n; ++i) {
            y[i*incY] *= x[i*incX];
        }
    }
}

template <typename IndexType, typename Alpha, typename VX, typename VY,
          typename Beta, typename VZ>
void
hvp(IndexType      n,
    const Alpha    &alpha,
    const VX       *x,
    IndexType      incX,
    const VY       *y,
    IndexType      incY,
    const Beta     &beta,
    VZ             *z,
    IndexType      incZ)
{
    const IndexType  UnitStride(1);

    if (incX==UnitStride && incY==UnitStride && incZ==UnitStride) {
        for (IndexType i=0; i<n; ++i) {
            z[i] = beta*z[i] + alpha*x[i]*y[i];
        }
    } else {
        for (IndexType i=0; i<n; ++i) {
            z[i*incZ] = beta*z[i*incZ] + alpha*x[i*incX]*y[i*incY];
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_HVP_TCC 1
