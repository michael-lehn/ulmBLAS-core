#ifndef ULMBLAS_CXXBLAS_LEVEL1_SWAP_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_SWAP_TCC 1

#include <ulmblas/cxxblas/level1/swap.h>
#include <ulmblas/impl/level1/swap.h>
#include <ulmblas/impl/level1extensions/geswap.h>

namespace cxxblas {

template <typename IndexType, typename VX, typename VY>
void
swap(IndexType      n,
     VX             *x,
     IndexType      incX,
     VY             *y,
     IndexType      incY)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }
    ulmBLAS::swap(n, x, incX, y, incY);
}

template <typename IndexType, typename MX, typename MY>
void
geswap(IndexType      m,
       IndexType      n,
       const MX       *X,
       IndexType      incRowX,
       IndexType      incColX,
       MY             *Y,
       IndexType      incRowY,
       IndexType      incColY)
{
    ulmBLAS::geswap(m, n, X, incRowX, incColX, Y, incRowY, incColY);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_SWAP_TCC
