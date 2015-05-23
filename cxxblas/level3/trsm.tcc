#ifndef ULMBLAS_CXXBLAS_LEVEL3_TRSM_TCC
#define ULMBLAS_CXXBLAS_LEVEL3_TRSM_TCC 1

#include <ulmblas/cxxblas/level3/trmm.h>
#include <ulmblas/external.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB>
void
trsm(bool         leftA,
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
    if (leftA) {
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

#ifdef WITH_CBLAS

template <typename IndexType, typename T>
typename cblas::If<IndexType,T>::isBlasCompatible
trsm(bool         leftA,
     IndexType    m,
     IndexType    n,
     const T      &alpha,
     bool         lowerA,
     bool         transA,
     bool         conjA,
     bool         unitDiagA,
     const T      *A,
     IndexType    incRowA,
     IndexType    incColA,
     T            *B,
     IndexType    incRowB,
     IndexType    incColB)
{
    const bool isColMajor = (incRowA==1) && (incRowB==1);
    const bool isRowMajor = (incColA==1) && (incColB==1);

    if (isColMajor || isRowMajor) {
        cblas::trsm(isColMajor, leftA, lowerA, transA, conjA, unitDiagA,
                    m, n,
                    alpha,
                    A, (isColMajor) ? incColA : incRowA,
                    B, (isColMajor) ? incColB : incRowB);
    } else {
        // fall back to ulmBLAS (or BLIS in future)
        trsm<IndexType, T, T, T>(leftA, m, n,
                                 alpha,
                                 lowerA, transA, conjA, unitDiagA,
                                 A, incRowA, incColA,
                                 B, incRowB, incColB);
    }
}

#endif // WITH_CBLAS

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_TRSM_TCC
