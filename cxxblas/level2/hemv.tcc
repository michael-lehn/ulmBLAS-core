#ifndef ULMBLAS_CXXBLAS_LEVEL2_HEMV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_HEMV_TCC 1

#include <ulmblas/cxxblas/level2/hemv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
hemv(IndexType    n,
     const Alpha  &alpha,
     bool         lowerA,
     const TA     *A,
     IndexType    incRowA,
     IndexType    incColA,
     const TX     *x,
     IndexType    incX,
     const Beta   &beta,
     TY           *y,
     IndexType    incY)
{
    if (incX<0) {
        x -= incX*(n-1);
    }
    if (incY<0) {
        y -= incY*(n-1);
    }

    if (lowerA) {
        ulmBLAS::helmv(n, alpha, false,
                        A, incRowA, incColA,
                        x, incX,
                        beta,
                        y, incY);
    } else {
        ulmBLAS::helmv(n, alpha, true,
                       A, incColA, incRowA,
                       x, incX,
                       beta,
                       y, incY);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_HEMV_TCC
