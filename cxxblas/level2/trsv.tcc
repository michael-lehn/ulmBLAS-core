#ifndef ULMBLAS_CXXBLAS_LEVEL2_TRSV_TCC
#define ULMBLAS_CXXBLAS_LEVEL2_TRSV_TCC 1

#include <ulmblas/cxxblas/level2/trmv.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
void
trsv(IndexType    n,
     bool         lowerA,
     bool         transA,
     bool         conjA,
     bool         unitDiagA,
     const TA     *A,
     IndexType    incRowA,
     IndexType    incColA,
     TX           *x,
     IndexType    incX)
{
    if (incX<0) {
        x -= incX*(n-1);
    }

//
//  Start the operations.
//
    if (lowerA) {
        if (!transA) {
            ulmBLAS::trlsv(n, unitDiagA, conjA, A, incRowA, incColA, x, incX);
        } else {
            ulmBLAS::trusv(n, unitDiagA, conjA, A, incColA, incRowA, x, incX);
        }
    } else {
        if (!transA) {
            ulmBLAS::trusv(n, unitDiagA, conjA, A, incRowA, incColA, x, incX);
        } else {
            ulmBLAS::trlsv(n, unitDiagA, conjA, A, incColA, incRowA, x, incX);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL2_TRSV_TCC
