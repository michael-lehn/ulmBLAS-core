#ifndef ULMBLAS_CXXBLAS_LEVEL2_SYR2_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_SYR2_TCC 1

#include <ulmblas/cxxblas/level2/her.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
void
syr2(IndexType    n,
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
        ulmBLAS::sylr2(n, alpha, x, incX, y, incY, A, incRowA, incColA);
    } else {
        ulmBLAS::sylr2(n, alpha, x, incX, y, incY, A, incColA, incRowA);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_SYR2_TCC
