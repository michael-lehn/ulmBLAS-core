#ifndef ULMBLAS_CXXBLAS_LEVEL1_CONJ_H
#define ULMBLAS_CXXBLAS_LEVEL1_CONJ_H 1

#include <ulmblas/cxxblas/level1/conj.h>
#include <ulmblas/impl/level1extensions/conj.h>

namespace cxxblas {

template <typename IndexType, typename VX>
void
conj(IndexType      n,
     VX             *x,
     IndexType      incX)
{
    ulmBLAS::conj(n, x, incX);
}

} // namespace cxxblas

#endif // ULMBLAS_CXXBLAS_LEVEL1_CONJ_H 1
