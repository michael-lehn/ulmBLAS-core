#ifndef ULMBLAS_CXXBLAS_LEVEL3_GEMM_TCC
#define ULMBLAS_CXXBLAS_LEVEL3_GEMM_TCC 1

#include <ulmblas/cxxblas/level3/gemm.h>
#include <ulmblas/external.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB,
          typename Beta, typename TC>
void
gemm(IndexType    m,
     IndexType    n,
     IndexType    k,
     const Alpha  &alpha,
     bool         transA,
     bool         conjA,
     const TA     *A,
     IndexType    incRowA_,
     IndexType    incColA_,
     bool         transB,
     bool         conjB,
     const TB     *B,
     IndexType    incRowB_,
     IndexType    incColB_,
     const Beta   &beta,
     TC           *C,
     IndexType    incRowC,
     IndexType    incColC)
{
    const IndexType incRowA = (!transA) ? incRowA_ : incColA_;
    const IndexType incColA = (!transA) ? incColA_ : incRowA_;

    const IndexType incRowB = (!transB) ? incRowB_ : incColB_;
    const IndexType incColB = (!transB) ? incColB_ : incRowB_;

    ulmBLAS::gemm(m, n, k,
                  alpha,
                  conjA, A, incRowA, incColA,
                  conjB, B, incRowB, incColB,
                  beta,
                  C, incRowC, incColC);
}

#ifdef WITH_CBLAS

template <typename IndexType, typename T>
typename cblas::If<IndexType,T>::isBlasCompatible
gemm(IndexType    m,
     IndexType    n,
     IndexType    k,
     const T      &alpha,
     bool         transA,
     bool         conjA,
     const T      *A,
     IndexType    incRowA,
     IndexType    incColA,
     bool         transB,
     bool         conjB,
     const T      *B,
     IndexType    incRowB,
     IndexType    incColB,
     const T      &beta,
     T            *C,
     IndexType    incRowC,
     IndexType    incColC)
{
    const bool isColMajor = (incRowA==1) && (incRowB==1) && (incRowC==1);
    const bool isRowMajor = (incColA==1) && (incColB==1) && (incColC==1);

    if (isColMajor || isRowMajor) {
        cblas::gemm(isColMajor, transA, conjA, transB, conjA,
                    m, n, k,
                    alpha,
                    A, (isColMajor) ? incColA : incRowA,
                    B, (isColMajor) ? incColB : incRowB,
                    beta,
                    C, (isColMajor) ? incColC : incRowC);
    } else {
        // fall back to ulmBLAS (or BLIS in future)
        gemm<IndexType, T, T, T, T, T>(m, n, k,
                                       alpha,
                                       transA, conjA, A, incRowA, incColA,
                                       transB, conjB, B, incRowB, incColB,
                                       beta,
                                       C, incRowC, incColC);
    }
}

#endif // WITH_CBLAS

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_GEMM_TCC
