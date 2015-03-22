#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEHMP_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEHMP_H 1

namespace ulmBLAS {

//
//  In-place Hadamard matrix product B = alpha * A .* B (component wise product)
//
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
           IndexType      incColB);

//
//  Hadamard matrix product C = beta*C + alpha * A .* B (component wise product)
//
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
          IndexType      incColC);

} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/gehmp.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GEHMP_H 1
