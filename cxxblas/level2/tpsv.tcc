#ifndef ULMBLAS_CXXBLAS_LEVEL2_TPSV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_TPSV_TCC 1

#include <ulmblas/cxxblas/level2/tpsv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
void
tpsv(IndexType    n,
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
            ulmBLAS::tplsv(n, unitDiagA, conjA, AP, x, incX);
        } else {
            ulmBLAS::tpusv(n, unitDiagA, conjA, AP, x, incX);
        }
    } else {
        if (lowerA) {
            ulmBLAS::tplstv(n, unitDiagA, conjA, AP, x, incX);
        } else {
            ulmBLAS::tpustv(n, unitDiagA, conjA, AP, x, incX);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_TPSV_TCC
