#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GERAXPY_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GERAXPY_H 1

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename MX, typename MY>
    void
    geraxpy(IndexType      m,
            IndexType      n,
            const Alpha    &alpha,
            const MX       *X,
            IndexType      incRowX,
            IndexType      incColX,
            MY             *Y,
            IndexType      incRowY,
            IndexType      incColY);

template <typename IndexType, typename Alpha, typename MX, typename MY>
    void
    geracxpy(IndexType      m,
             IndexType      n,
             const Alpha    &alpha,
             const MX       *X,
             IndexType      incRowX,
             IndexType      incColX,
             MY             *Y,
             IndexType      incRowY,
             IndexType      incColY);

} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/geaxpy.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GERAXPY_H 1
