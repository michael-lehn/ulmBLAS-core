#ifndef ULMBLAS_CXXBLAS_LEVEL1_IAMAX_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_IAMAX_TCC 1

#include <ulmblas/cxxblas/level1/iamax.h>
#include <ulmblas/impl/level1/iamax.h>

namespace cxxblas {

template <typename IndexType, typename VX>
IndexType
iamax(IndexType      n,
      const VX       *x,
      IndexType      incX)
{
    return ulmBLAS::iamax(n, x, incX);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_IAMAX_TCC
