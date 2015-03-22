#ifndef ULMBLAS_CXXBLAS_LEVEL2_TPMV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_TPMV_TCC 1

#include <ulmblas/cxxblas/level2/tpmv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
void
tpmv(IndexType    n,
     bool         colMajorA,
     bool         lowerA,
     bool         transA,
     bool         conjA,
     bool         unitDiagA,
     const TA     *AP,
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
            ulmBLAS::tplmv(n, unitDiagA, conjA, AP, x, incX);
        } else {
            ulmBLAS::tpumv(n, unitDiagA, conjA, AP, x, incX);
        }
    } else {
        if (lowerA) {
            ulmBLAS::tplmtv(n, unitDiagA, conjA, AP, x, incX);
        } else {
            ulmBLAS::tpumtv(n, unitDiagA, conjA, AP, x, incX);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_TPMV_TCC
