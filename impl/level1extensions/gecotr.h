#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_H 1

namespace ulmBLAS {

template <typename IndexType, typename MX>
    void
    gecotr(IndexType      n,
           bool           transX,
           bool           conjX,
           MX             *X,
           IndexType      incRowX,
           IndexType      incColX);


} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/gecotr.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_H 1
