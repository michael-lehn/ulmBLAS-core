#ifndef ULMBLAS_CXXBLAS_LEVEL1_HVP_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_HVP_TCC 1

#include <ulmblas/cxxblas/level1/hvp.h>
#include <ulmblas/impl/level1extensions/hvp.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename VX, typename VY>
void
ihvp(IndexType      n,
     const Alpha    &alpha,
     const VX       *x,
     IndexType      incX,
     VY             *y,
     IndexType      incY)
{
    ulmBLAS::ihvp(n, alpha, x, incX, y, incY);
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
    ulmBLAS::hvp(n, alpha, x, incX, y, incY, beta, z, incZ);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_HVP_TCC
