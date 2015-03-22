#ifndef ULMBLAS_CXXBLAS_LEVEL2_TBMV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_TBMV_TCC 1

#include <ulmblas/cxxblas/level2/tbmv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
void
tbmv(IndexType    n,
     IndexType    k,
     bool         colMajorA,
     bool         lowerA,
     bool         transA,
     bool         conjA,
     bool         unitDiagA,
     const TA     *A,
     IndexType    ldA,
     TX           *x,
     IndexType    incX)
{
    if (!colMajorA) {
        lowerA = !lowerA;
        transA = !transA;
    }

    if (incX<0) {
        x -= incX*(n-1);
    }

    if (!transA) {
        if (lowerA) {
            ulmBLAS::tblmv(n, k, unitDiagA, conjA, A, ldA, x, incX);
        } else {
            ulmBLAS::tbumv(n, k, unitDiagA, conjA, A, ldA, x, incX);
        }
    } else {
        if (lowerA) {
            ulmBLAS::tblmtv(n, k, unitDiagA, conjA, A, ldA, x, incX);
        } else {
            ulmBLAS::tbumtv(n, k, unitDiagA, conjA, A, ldA, x, incX);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_TBMV_TCC
