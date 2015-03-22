#ifndef ULMBLAS_CXXBLAS_LEVEL1_GEHMP_H
#define ULMBLAS_CXXBLAS_LEVEL1_GEHMP_H 1

namespace cxxblas {

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

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/gehmp.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_GEHMP_H 1
