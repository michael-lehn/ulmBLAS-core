#ifndef ULMBLAS_CXXBLAS_LEVEL2_TBSV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_TBSV_TCC 1

#include <ulmblas/cxxblas/level2/tbsv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
void
tbsv(IndexType    n,
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
            ulmBLAS::tblsv(n, k, unitDiagA, conjA, A, ldA, x, incX);
        } else {
            ulmBLAS::tbusv(n, k, unitDiagA, conjA, A, ldA, x, incX);
        }
    } else {
        if (lowerA) {
            ulmBLAS::tblstv(n, k, unitDiagA, conjA, A, ldA, x, incX);
        } else {
            ulmBLAS::tbustv(n, k, unitDiagA, conjA, A, ldA, x, incX);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_TBSV_TCC
