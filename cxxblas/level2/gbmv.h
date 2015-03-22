#ifndef ULMBLAS_CXXBLAS_LEVEL2_GBMV_H
#define ULMBLAS_CXXBLAS_LEVEL2_GBMV_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TA, typename TX,
          typename Beta, typename TY>
    void
    gbmv(IndexType    m,
         IndexType    n,
         IndexType    kl,
         IndexType    ku,
         const Alpha  &alpha,
         bool         colMajorA,
         bool         transA,
         bool         conjA,
         const TA     *A,
         IndexType    ldA,
         const TX     *x,
         IndexType    incX,
         const Beta   &beta,
         TY           *y,
         IndexType    incY);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/gbmv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_GBMV_H
