#ifndef ULMBLAS_CXXBLAS_LEVEL1_GECOTR_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_GECOTR_TCC 1

#include <ulmblas/cxxblas/level1/gecotr.h>

namespace cxxblas {

template <typename IndexType, typename MX, typename MY>
void
gecotr(IndexType      n,
       bool           transX,
       bool           conjX,
       MX             *X,
       IndexType      incRowX,
       IndexType      incColX)
{
    ulmBLAS::gecotr(n, transX, conjX, X, incRowX, incColX);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_GECOTR_TCC
