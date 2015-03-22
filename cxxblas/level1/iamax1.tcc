#ifndef ULMBLAS_CXXBLAS_LEVEL1_IAMAX1_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_IAMAX1_TCC 1

#include <ulmblas/cxxblas/level1/iamax1.h>
#include <ulmblas/impl/level1extensions/iamax1.h>

namespace cxxblas {

template <typename IndexType, typename VX>
IndexType
iamax1(IndexType      n,
       const VX       *x,
       IndexType      incX)
{
    return ulmBLAS::iamax1(n, x, incX);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_IAMAX1_TCC
