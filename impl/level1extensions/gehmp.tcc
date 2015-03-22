#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEHMP_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEHMP_TCC 1

#include <ulmblas/impl/level1extensions/gehmp.h>
#include <ulmblas/impl/level1extensions/hvp.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename MA, typename MB>
void
geihmp(IndexType      m,
       IndexType      n,
       const Alpha    &alpha,
       const MA       *A,
       IndexType      incRowA,
       IndexType      incColA,
       MB             *B,
       IndexType      incRowB,
       IndexType      incColB)
{
    if (incRowA<incColA) {
        for (IndexType j=0; j<n; ++j) {
            hvp(m, alpha, &A[j*incColA], incRowA, &B[j*incColB], incRowB);
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            hvp(n, alpha, &A[i*incRowA], incColA, &B[i*incRowB], incColB);
        }
    }
}

template <typename IndexType, typename Alpha, typename MA, typename MB,
          typename Beta, typename MC>
void
gehmp(IndexType      m,
      IndexType      n,
      const Alpha    &alpha,
      const MA       *A,
      IndexType      incRowA,
      IndexType      incColA,
      const MB       *B,
      IndexType      incRowB,
      IndexType      incColB,
      const Beta     &beta,
      MC             *C,
      IndexType      incRowC,
      IndexType      incColC)
{
    if (incRowA<incColA) {
        for (IndexType j=0; j<n; ++j) {
            hvp(m,
                alpha,
                &A[j*incColA], incRowA,
                &B[j*incColB], incRowB,
                beta,
                &C[j*incColC], incRowC);
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            hvp(n,
                alpha,
                &A[i*incRowA], incColA,
                &B[i*incRowB], incColB,
                beta,
                &C[i*incRowC], incColC);
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEHMP_TCC 1
