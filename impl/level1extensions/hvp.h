#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_HVP_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_HVP_H 1

namespace ulmBLAS {

//
//  In-place Hadamard vector product y = alpha * x .* y (component wise product)
//
template <typename IndexType, typename Alpha, typename VX, typename VY>
    void
    ihvp(IndexType      n,
         const Alpha    &alpha,
         const VX       *x,
         IndexType      incX,
         VY             *y,
         IndexType      incY);

//
//  Hadamard vector product z = beta*z + alpha * x .* y (component wise product)
//
template <typename IndexType, typename Alpha, typename VX, typename VY,
          typename Beta, typename VZ>
    void
    hvp(IndexType      n,
        const Alpha    &alpha,
        const VX       *x,
        IndexType      incX,
        const VY       *y,
        IndexType      incY,
        const Beta     &beta,
        VZ             *z,
        IndexType      incZ);

} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/hvp.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_HVP_H 1
