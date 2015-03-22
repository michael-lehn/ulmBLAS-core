#ifndef ULMBLAS_CXXBLAS_LEVEL2_GEMV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_GEMV_TCC 1

#include <ulmblas/cxxblas/level2/gemv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
gemv(IndexType    m,
     IndexType    n,
     const Alpha  &alpha,
     bool         transA,
     bool         conjA,
     const TA     *A,
     IndexType    incRowA,
     IndexType    incColA,
     const TX     *x,
     IndexType    incX,
     const Beta   &beta,
     TY           *y,
     IndexType    incY)
{
    if (!transA) {
        if (incX<0) {
            x -= incX*(n-1);
        }
        if (incY<0) {
            y -= incY*(m-1);
        }
        ulmBLAS::gemv(m, n,
                      alpha,
                      conjA, A, incRowA, incColA,
                      x, incX,
                      beta,
                      y, incY);
    } else {
        if (incX<0) {
            x -= incX*(m-1);
        }
        if (incY<0) {
            y -= incY*(n-1);
        }
        ulmBLAS::gemv(n, m,
                      alpha,
                      conjA, A, incColA, incRowA,
                      x, incX,
                      beta,
                      y, incY);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_GEMV_TCC
