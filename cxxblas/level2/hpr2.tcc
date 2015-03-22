#ifndef ULMBLAS_CXXBLAS_LEVEL2_HPR2_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_HPR2_TCC 1

#include <ulmblas/cxxblas/level2/hpr2.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
void
hpr2(IndexType    n,
     const Alpha  &alpha,
     const TX     *x,
     IndexType    incX,
     const TY     *y,
     IndexType    incY,
     bool         colMajorA,
     bool         lowerA,
     TA           *A)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }

    if (colMajorA) {
        if (lowerA) {
            ulmBLAS::hplr2(n, alpha, x, incX, y, incY, A);
        } else {
            ulmBLAS::hpur2(n, alpha, x, incX, y, incY, A);
        }
    } else {
        if (!lowerA) {
            ulmBLAS::hplr2(n, true, alpha, x, incX, y, incY, A);
        } else {
            ulmBLAS::hpur2(n, true, alpha, x, incX, y, incY, A);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_HPR2_TCC
