#ifndef ULMBLAS_CXXBLAS_LEVEL2_GERC_H
#define ULMBLAS_CXXBLAS_LEVEL2_GERC_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
    void
    gerc(IndexType    m,
         IndexType    n,
         const Alpha  &alpha,
         const TX     *x,
         IndexType    incX,
         const TY     *y,
         IndexType    incY,
         TA           *A,
         IndexType    incRowA,
         IndexType    incColA);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/gerc.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_GERC_H

