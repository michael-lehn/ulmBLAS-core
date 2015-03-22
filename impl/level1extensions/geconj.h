#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECONJ_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECONJ_H 1

namespace ulmBLAS {

template <typename IndexType, typename MX>
    void
    geconj(IndexType      m,
           IndexType      n,
           MX             *X,
           IndexType      incRowX,
           IndexType      incColX);

} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/geconj.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECONJ_H 1
