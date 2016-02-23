#ifndef ULMBLAS_CXXBLAS_LEVEL2_GERU_H
#define ULMBLAS_CXXBLAS_LEVEL2_GERU_H 1

#ifdef USE_EXTERNAL_BLAS
#include <ulmblas/external/blis/geru.h>
#endif // USE_EXTERNAL_BLAS

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TY,
          typename TA>
    void
    geru(IndexType    m,
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

#include <ulmblas/cxxblas/level2/geru.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_GERU_H

