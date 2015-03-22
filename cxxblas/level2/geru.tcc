#ifndef ULMBLAS_CXXBLAS_LEVEL2_GERU_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_GERU_TCC 1

#include <ulmblas/cxxblas/level2/geru.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
void
geru(IndexType    m,
     IndexType    n,
     const Alpha  &alpha,
     const TX     *x,
     IndexType    incX,
     const TY     *y,
     IndexType    incY,
     TA           *A,
     IndexType    incRowA,
     IndexType    incColA)
{
    if (incX<0) {
        x -= incX*(m-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }

    ulmBLAS::ger(m, n, alpha, x, incX, y, incY, A, incRowA, incColA);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_GERU_TCC
