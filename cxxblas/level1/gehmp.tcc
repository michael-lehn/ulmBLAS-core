#ifndef ULMBLAS_CXXBLAS_LEVEL1_GEHMP_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_GEHMP_TCC 1

#include <ulmblas/cxxblas/level1/gehmp.h>
#include <ulmblas/impl/level1extensions/gehmp.h>

namespace cxxblas {

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
    ulmBLAS::geihmp(m, n, alpha, A, incRowA, incColA, B, incRowB, incColB);
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
    ulmBLAS::gehmp(m, n,
                   alpha,
                   A, incRowA, incColA,
                   B, incRowB, incColB,
                   beta,
                   C, incRowC, incColC);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_GEHMP_TCC 1
