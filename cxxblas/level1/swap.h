#ifndef ULMBLAS_CXXBLAS_LEVEL1_SWAP_H
#define ULMBLAS_CXXBLAS_LEVEL1_SWAP_H 1

namespace cxxblas {

template <typename IndexType, typename VX, typename VY>
    void
    swap(IndexType      n,
         VX             *x,
         IndexType      incX,
         VY             *y,
         IndexType      incY);

template <typename IndexType, typename MX, typename MY>
    void
    geswap(IndexType      m,
           IndexType      n,
           const MX       *X,
           IndexType      incRowX,
           IndexType      incColX,
           MY             *Y,
           IndexType      incRowY,
           IndexType      incColY);

} // namespace cxxblas

#include <ulmblas/cxxblas/level1/swap.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL1_SWAP_H
