#ifndef ULMBLAS_CXXBLAS_LEVEL3_SYR2K_TCC
#define ULMBLAS_CXXBLAS_LEVEL3_SYR2K_TCC 1

#include <ulmblas/cxxblas/level3/syr2k.h>
#include <ulmblas/impl/level3/sylr2k.h>
#include <ulmblas/impl/level3/syur2k.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TB,
          typename Beta, typename TC>
void
syr2k(IndexType    n,
      IndexType    k,
      const Alpha  &alpha,
      bool         transAB,
      const TA     *A,
      IndexType    incRowA,
      IndexType    incColA,
      const TB     *B,
      IndexType    incRowB,
      IndexType    incColB,
      const Beta   &beta,
      bool         lowerC,
      TC           *C,
      IndexType    incRowC,
      IndexType    incColC)
{
    if (!transAB) {
        if (lowerC) {
            ulmBLAS::sylr2k(n, k,
                            alpha,
                            A, incRowA, incColA,
                            B, incRowB, incColB,
                            beta,
                            C, incRowC, incColC);
        } else {
            ulmBLAS::syur2k(n, k,
                            alpha,
                            A, incRowA, incColA,
                            B, incRowB, incColB,
                            beta,
                            C, incRowC, incColC);
        }
    } else {
        if (lowerC) {
            ulmBLAS::syur2k(n, k,
                            alpha,
                            A, incColA, incRowA,
                            B, incColB, incRowB,
                            beta,
                            C, incColC, incRowC);
        } else {
            ulmBLAS::sylr2k(n, k,
                            alpha,
                            A, incColA, incRowA,
                            B, incColB, incRowB,
                            beta,
                            C, incColC, incRowC);
        }
    }
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL3_SYR2K_TCC
