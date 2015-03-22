#ifndef ULMBLAS_CXXBLAS_LEVEL1_RSCAL_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_RSCAL_TCC 1

#include <ulmblas/cxxblas/level1/rscal.h>
#include <ulmblas/impl/level1extensions/rscal.h>
#include <ulmblas/impl/level1extensions/gerscal.h>

namespace cxxblas {

template <typename IndexType, typename Alpha, typename VX>
void
rscal(IndexType      n,
      const Alpha    &alpha,
      VX             *x,
      IndexType      incX)
{
    ulmBLAS::rscal(n, alpha, x, incX);
}

template <typename IndexType, typename Alpha, typename MA>
void
gerscal(IndexType    m,
        IndexType    n,
        const Alpha  &alpha,
        MA           *A,
        IndexType    incRowA,
        IndexType    incColA)
{
    ulmBLAS::gerscal(m, n, alpha, A, incRowA, incColA);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_RSCAL_TCC
