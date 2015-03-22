#ifndef ULMBLAS_CXXBLAS_LEVEL1_ROT_H
#define ULMBLAS_CXXBLAS_LEVEL1_ROT_H 1

namespace cxxblas {

template <typename A, typename B, typename T>
    void
    rotg(A &a,
         B &b,
         T &c,
         T &s);

template <typename TA, typename TB, typename T>
    void
    rotg(std::complex<TA>   &a,
         std::complex<TB>   &b,
         T                  &c,
         std::complex<T>    &s);

template <typename IndexType, typename VX, typename VY, typename T>
    void
    rot(IndexType   n,
        VX          *x,
        IndexType   incX,
        VY          *y,
        IndexType   incY,
        T           c,
        T           s);

template <typename IndexType, typename X, typename Y, typename T>
    void
    rot(IndexType              n,
        std::complex<X>        *x,
        IndexType              incX,
        std::complex<Y>        *y,
        IndexType              incY,
        T                      c,
        const std::complex<T>  &s);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/rot.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_ROT_H
