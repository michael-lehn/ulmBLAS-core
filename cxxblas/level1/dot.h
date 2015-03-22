#ifndef ULMBLAS_CXXBLAS_LEVEL1_DOT_H
#define ULMBLAS_CXXBLAS_LEVEL1_DOT_H 1

namespace cxxblas {

template <typename IndexType, typename VX, typename VY, typename Result>
    void
    dotu(IndexType      n,
         const VX       *x,
         IndexType      incX,
         const VY       *y,
         IndexType      incY,
         Result         &result);

template <typename IndexType, typename VX, typename VY, typename Result>
    void
    dotc(IndexType      n,
         const VX       *x,
         IndexType      incX,
         const VY       *y,
         IndexType      incY,
         Result         &result);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/dot.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_DOT_H
