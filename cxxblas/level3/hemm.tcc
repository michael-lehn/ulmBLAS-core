#ifndef ULMBLAS_CXXBLAS_LEVEL3_HEMM_TCC
#define ULMBLAS_CXXBLAS_LEVEL3_HEMM_TCC 1

#include <ulmblas/cxxblas/level3/hemm.tcc>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB,
          typename Beta, typename TC>
void
hemm(bool         left,
     IndexType    m,
     IndexType    n,
     const Alpha  &alpha,
     bool         lowerA,
     const TA     *A,
     IndexType    incRowA,
     IndexType    incColA,
     const TB     *B,
     IndexType    incRowB,
     IndexType    incColB,
     const Beta   &beta,
     TC           *C,
     IndexType    incRowC,
     IndexType    incColC)
{
    if (left) {
        if (lowerA) {
            ulmBLAS::helmm(m, n,
                           alpha,
                           A, incRowA, incColA,
                           B, incRowB, incColB,
                           beta,
                           C, incRowC, incColC);
        } else {
            ulmBLAS::heumm(m, n,
                           alpha,
                           A, incRowA, incColA,
                           B, incRowB, incColB,
                           beta,
                           C, incRowC, incColC);
        }
    } else {
        if (lowerA) {
            ulmBLAS::heumm(n, m,
                           alpha,
                           A, incColA, incRowA,
                           B, incColB, incRowB,
                           beta,
                           C, incColC, incRowC);
        } else {
            ulmBLAS::helmm(n, m,
                           alpha,
                           A, incColA, incRowA,
                           B, incColB, incRowB,
                           beta,
                           C, incColC, incRowC);
        }
    }

}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_HEMM_TCC
