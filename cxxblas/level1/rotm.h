#ifndef ULMBLAS_CXXBLAS_LEVEL1_ROTM_H
#define ULMBLAS_CXXBLAS_LEVEL1_ROTM_H 1

namespace cxxblas {

template <typename IndexType, typename T>
    void
    rotmg(T &d1, T &d2, T &x1, T &y1, T *param);

template <typename IndexType, typename T>
    void
    rotm(IndexType  n,
         T          *x,
         IndexType  incX,
         T          *y,
         IndexType  incY,
         T          *param);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/rot.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_ROTM_H
