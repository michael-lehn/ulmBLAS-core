#ifndef ULMBLAS_CXXBLAS_LEVEL1_NRM2_TCC
#define ULMBLAS_CXXBLAS_LEVEL1_NRM2_TCC 1

#include <ulmblas/cxxblas/level1/nrm2.h>
#include <ulmblas/impl/level1/nrm2.h>

namespace cxxblas {

template <typename IndexType, typename VX, typename Result>
void
nrm2(IndexType  n,
     const VX   *x,
     IndexType  incX,
     Result     &result)
{
    ulmBLAS::nrm2(n, x, incX, result);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_NRM2_TCC
