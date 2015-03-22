#ifndef ULMBLAS_CXXBLAS_LEVEL1_DOT_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_DOT_TCC 1

#include <ulmblas/cxxblas/level1/dot.h>
#include <ulmblas/impl/level1/dot.h>

namespace cxxblas {

template <typename IndexType, typename VX, typename VY, typename Result>
void
dotu(IndexType      n,
     const VX       *x,
     IndexType      incX,
     const VY       *y,
     IndexType      incY,
     Result         &result)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }
    ulmBLAS::dotu(n, x, incX, y, incY, result);
}

template <typename IndexType, typename VX, typename VY, typename Result>
void
dotc(IndexType      n,
     const VX       *x,
     IndexType      incX,
     const VY       *y,
     IndexType      incY,
     Result         &result)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }
    ulmBLAS::dotc(n, x, incX, y, incY, result);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_DOT_TCC
