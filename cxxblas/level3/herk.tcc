#ifndef ULMBLAS_CXXBLAS_LEVEL3_HERK_TCC
#define ULMBLAS_CXXBLAS_LEVEL3_HERK_TCC 1

#include <ulmblas/cxxblas/level3/herk.h>
#include <ulmblas/ulmblas.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename Beta,
          typename TC>
void
herk(IndexType    n,
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
            ulmBLAS::helrk(n, k,
                           alpha,
                           A, incRowA, incColA,
                           beta,
                           C, incRowC, incColC);
        } else {
            ulmBLAS::heurk(n, k,
                           alpha,
                           A, incRowA, incColA,
                           beta,
                           C, incRowC, incColC);
        }
    } else {
        if (lowerC) {
            ulmBLAS::heurk(n, k,
                           alpha,
                           A, incColA, incRowA,
                           beta,
                           C, incColC, incRowC);
        } else {
            ulmBLAS::helrk(n, k,
                           alpha,
                           A, incColA, incRowA,
                           beta,
                           C, incColC, incRowC);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_HERK_TCC

