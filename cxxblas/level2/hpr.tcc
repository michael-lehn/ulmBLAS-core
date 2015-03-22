#ifndef ULMBLAS_CXXBLAS_LEVEL2_HPR_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_HPR_TCC 1

#include <ulmblas/cxxblas/level2/hpr.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TA>
void
hpr(IndexType    n,
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
            ulmBLAS::hplr(n, alpha, x, incX, A);
        } else {
            ulmBLAS::hpur(n, alpha, x, incX, A);
        }
    } else {
        if (!lowerA) {
            ulmBLAS::hplr(n, alpha, true, x, incX, A);
        } else {
            ulmBLAS::hpur(n, alpha, true, x, incX, A);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_HPR_TCC
