#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_H 1

namespace ulmBLAS {

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

template <typename IndexType, typename MX, typename MY>
    void
    geswapc(IndexType      m,
            IndexType      n,
            const MX       *X,
            IndexType      incRowX,
            IndexType      incColX,
            MY             *Y,
            IndexType      incRowY,
            IndexType      incColY);


} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/geswap.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_H 1
