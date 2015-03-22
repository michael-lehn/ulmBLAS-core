#ifndef ULMBLAS_CXXBLAS_LEVEL2_SYR_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_SYR_TCC 1

#include <ulmblas/cxxblas/level2/her.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TA>
void
syr(IndexType    n,
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
        ulmBLAS::sylr(n, alpha, x, incX, A, incRowA, incColA);
    } else {
        ulmBLAS::sylr(n, alpha, x, incX, A, incColA, incRowA);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_SYR_TCC
