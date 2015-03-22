#ifndef ULMBLAS_CXXBLAS_LEVEL2_HER2_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_HER2_TCC 1

#include <ulmblas/cxxblas/level2/her.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
void
her2(IndexType    n,
     const Alpha  &alpha,
     const TX     *x,
     IndexType    incX,
     const TY     *y,
     IndexType    incY,
     bool         lowerA,
     TA           *A,
     IndexType    incRowA,
     IndexType    incColA)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }

    if (lowerA) {
        ulmBLAS::helr2(n, false, alpha, x, incX, y, incY, A, incRowA, incColA);
    } else {
        ulmBLAS::helr2(n, true, alpha, x, incX, y, incY, A, incColA, incRowA);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_HER2_TCC
