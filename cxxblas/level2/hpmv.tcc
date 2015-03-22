#ifndef ULMBLAS_CXXBLAS_LEVEL2_HPMV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_HPMV_TCC 1

#include <ulmblas/cxxblas/level2/hpmv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
hpmv(IndexType    n,
     const Alpha  &alpha,
     bool         colMajorA,
     bool         lowerA,
     const TA     *AP,
     const TX     *x,
     IndexType    incX,
     const Beta   &beta,
     TY           *y,
     IndexType    incY)
{
    if (!colMajorA) {
        lowerA = !lowerA;
    }

    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }

    if (lowerA) {
        ulmBLAS::hplmv(n, alpha, AP, x, incX, beta, y, incY);
    } else {
        ulmBLAS::hpumv(n, alpha, AP, x, incX, beta, y, incY);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_HPMV_TCC
