#ifndef ULMBLAS_CXXBLAS_LEVEL2_HER_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_HER_TCC 1

#include <ulmblas/cxxblas/level2/her.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TA>
void
her(IndexType    n,
    const Alpha  &alpha,
    const TX     *x,
    IndexType    incX,
    bool         lowerA,
    TA           *A,
    IndexType    incRowA,
    IndexType    incColA)
{
    if (incX<0) {
        x -= incX*(n-1);
    }

    if (lowerA) {
        ulmBLAS::helr(n, alpha, false, x, incX, A, incRowA, incColA);
    } else {
        ulmBLAS::helr(n, alpha, true, x, incX, A, incColA, incRowA);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_HER_TCC
