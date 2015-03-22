#ifndef ULMBLAS_CXXBLAS_LEVEL3_TRSM_TCC
#define ULMBLAS_CXXBLAS_LEVEL3_TRSM_TCC 1

#include <ulmblas/cxxblas/level3/trmm.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB>
void
trsm(bool         left,
     IndexType    m,
     IndexType    n,
     const Alpha  &alpha,
     bool         lowerA,
     bool         transA,
     bool         conjA,
     bool         unitDiagA,
     const TA     *A,
     IndexType    incRowA,
     IndexType    incColA,
     TB           *B,
     IndexType    incRowB,
     IndexType    incColB)
{
    if (left) {
        if (lowerA) {
            if (!transA) {
                ulmBLAS::trlsm(m, n, alpha, conjA, unitDiagA,
                               A, incRowA, incColA,
                               B, incRowB, incColB);
            } else {
                ulmBLAS::trusm(m, n, alpha, conjA, unitDiagA,
                               A, incColA, incRowA,
                               B, incRowB, incColB);
            }
        } else {
            if (!transA) {
                ulmBLAS::trusm(m, n, alpha, conjA, unitDiagA,
                               A, incRowA, incColA,
                               B, incRowB, incColB);
            } else {
                ulmBLAS::trlsm(m, n, alpha, conjA, unitDiagA,
                               A, incColA, incRowA,
                               B, incRowB, incColB);
            }
        }
    } else {
        if (lowerA) {
            if (!transA) {
                ulmBLAS::trusm(n, m, alpha, conjA, unitDiagA,
                               A, incColA, incRowA,
                               B, incColB, incRowB);
            } else {
                ulmBLAS::trlsm(n, m, alpha, conjA, unitDiagA,
                               A, incRowA, incColA,
                               B, incColB, incRowB);
            }
        } else {
            if (!transA) {
                ulmBLAS::trlsm(n, m, alpha, conjA, unitDiagA,
                               A, incColA, incRowA,
                               B, incColB, incRowB);
            } else {
                ulmBLAS::trusm(n, m, alpha, conjA, unitDiagA,
                               A, incRowA, incColA,
                               B, incColB, incRowB);
            }
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_TRSM_TCC
