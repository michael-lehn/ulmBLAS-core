#ifndef ULMBLAS_CXXBLAS_LEVEL2_GEMV_H
#define ULMBLAS_CXXBLAS_LEVEL2_GEMV_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
    void
    gemv(IndexType    m,
         IndexType    n,
         const Alpha  &alpha,
         bool         transA,
         bool         conjA,
         const TA     *A,
         IndexType    incRowA,
         IndexType    incColA,
         const TX     *x,
         IndexType    incX,
         const Beta   &beta,
         TY           *y,
         IndexType    incY);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/gemv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_GEMV_H

