#ifndef ULMBLAS_CXXBLAS_LEVEL2_SPR_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_SPR_TCC 1

#include <ulmblas/cxxblas/level2/hpr.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TA>
void
spr(IndexType    n,
    const Alpha  &alpha,
    const TX     *x,
    IndexType    incX,
    bool         colMajorA,
    bool         lowerA,
    TA           *A)
{
    if (incX<0) {
        x -= incX*(n-1);
    }

    if (colMajorA) {
        if (lowerA) {
            ulmBLAS::splr(n, alpha, x, incX, A);
        } else {
            ulmBLAS::spur(n, alpha, x, incX, A);
        }
    } else {
        if (!lowerA) {
            ulmBLAS::splr(n, alpha, x, incX, A);
        } else {
            ulmBLAS::spur(n, alpha, x, incX, A);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_SPR_TCC
