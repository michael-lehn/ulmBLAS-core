#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPYF_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPYF_H 1

#include <ulmblas/impl/level1extensions/kernel/axpyf.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename VA, typename VX,
          typename VY>
    void
    axpyf(IndexType      n,
          const Alpha    &alpha,
          const VA       *a,
          IndexType      incA,
          const VX       *X,
          IndexType      incRowX,
          IndexType      incColX,
          VY             *y,
          IndexType      incY);

template <typename IndexType, typename Alpha, typename VA, typename VX,
          typename VY>
    void
    acxpyf(IndexType      n,
           const Alpha    &alpha,
           const VA       *a,
           IndexType      incA,
           const VX       *X,
           IndexType      incRowX,
           IndexType      incColX,
           VY             *y,
           IndexType      incY);

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_AXPYF_H

#include <ulmblas/impl/level1extensions/axpyf.tcc>
