#ifndef ULMBLAS_CXXBLAS_LEVEL2_SPR2_H
#define ULMBLAS_CXXBLAS_LEVEL2_SPR2_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
    void
    spr2(IndexType    n,
         const Alpha  &alpha,
         const TX     *x,
         IndexType    incX,
         const TY     *y,
         IndexType    incY,
         bool         colMajorA,
         bool         lowerA,
         TA           *A);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/spr2.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_SPR2_H
