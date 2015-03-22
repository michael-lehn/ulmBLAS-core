#ifndef ULMBLAS_CXXBLAS_LEVEL2_HBMV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_HBMV_TCC 1

#include <ulmblas/cxxblas/level2/hbmv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
hbmv(IndexType    n,
     IndexType    k,
     const Alpha  &alpha,
     bool         colMajorA,
     bool         lowerA,
     const TA     *A,
     IndexType    ldA,
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

//
//  Start the operations.
//
    if (lowerA) {
        ulmBLAS::hblmv(n, k, alpha, A, ldA, x, incX, beta, y, incY);
    } else {
        ulmBLAS::hbumv(n, k, alpha, A, ldA, x, incX, beta, y, incY);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_HBMV_TCC
