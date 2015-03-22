#ifndef ULMBLAS_CXXBLAS_LEVEL1_GECONJ_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_GECONJ_TCC 1

#include <ulmblas/cxxblas/level1/geconj.h>

namespace cxxblas {

template <typename IndexType, typename MX, typename MY>
void
geconj(IndexType      m,
       IndexType      n,
       MX             *X,
       IndexType      incRowX,
       IndexType      incColX)
{
    ulmBLAS::geconj(m, n, X, incRowX, incColX);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_GECONJ_TCC
