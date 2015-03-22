#ifndef ULMBLAS_CXXBLAS_LEVEL1_RSCAL_H
#define ULMBLAS_CXXBLAS_LEVEL1_RSCAL_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename VX>
    void
    rscal(IndexType      n,
          const Alpha    &alpha,
          VX             *x,
          IndexType      incX);

template <typename IndexType, typename Alpha, typename MA>
    void
    gerscal(IndexType    m,
            IndexType    n,
            const Alpha  &alpha,
            MA           *A,
            IndexType    incRowA,
            IndexType    incColA);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/rscal.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_RSCAL_H
