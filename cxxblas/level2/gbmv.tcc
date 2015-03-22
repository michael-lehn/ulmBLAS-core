#ifndef ULMBLAS_CXXBLAS_LEVEL2_GBMV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_GBMV_TCC 1

#include <ulmblas/cxxblas/level2/gbmv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
void
gbmv(IndexType    m,
     IndexType    n,
     IndexType    kl,
     IndexType    ku,
     const Alpha  &alpha,
     bool         colMajorA,
     bool         transA,
     bool         conjA,
     const TA     *A,
     IndexType    ldA,
     const TX     *x,
     IndexType    incX,
     const Beta   &beta,
     TY           *y,
     IndexType    incY)
{
    if (!colMajorA) {
        transA = !transA;
    }

    if (!transA) {
        if (incX<0) {
            x -= incX*(n-1);
        }
        if (incY<0) {
            y -= incY*(m-1);
        }
    } else {
        if (incX<0) {
            x -= incX*(m-1);
        }
        if (incY<0) {
            y -= incY*(n-1);
        }
    }

    if (!transA) {
        ulmBLAS::gbmv(m, n, kl, ku, alpha,
                      conjA, A, ldA,
                      x, incX,
                      beta,
                      y, incY);
    } else {
        ulmBLAS::gbmtv(m, n, kl, ku, alpha,
                       conjA, A, ldA,
                       x, incX,
                       beta,
                       y, incY);
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_GBMV_TCC
