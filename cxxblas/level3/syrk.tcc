#ifndef ULMBLAS_CXXBLAS_LEVEL3_SYRK_TCC
#define ULMBLAS_CXXBLAS_LEVEL3_SYRK_TCC 1

#include <ulmblas/cxxblas/level3/syrk.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename Beta,
          typename TC>
void
syrk(IndexType    n,
     IndexType    k,
     const Alpha  &alpha,
     bool         transA,
     const TA     *A,
     IndexType    incRowA,
     IndexType    incColA,
     const Beta   &beta,
     bool         lowerC,
     TC           *C,
     IndexType    incRowC,
     IndexType    incColC)
{
    if (!transA) {
        if (lowerC) {
            ulmBLAS::sylrk(n, k,
                           alpha,
                           A, incRowA, incColA,
                           beta,
                           C, incRowC, incColC);
        } else {
            ulmBLAS::syurk(n, k,
                           alpha,
                           A, incRowA, incColA,
                           beta,
                           C, incRowC, incColC);
        }
    } else {
        if (lowerC) {
            ulmBLAS::syurk(n, k,
                           alpha,
                           A, incColA, incRowA,
                           beta,
                           C, incColC, incRowC);
        } else {
            ulmBLAS::sylrk(n, k,
                           alpha,
                           A, incColA, incRowA,
                           beta,
                           C, incColC, incRowC);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_SYRK_TCC

