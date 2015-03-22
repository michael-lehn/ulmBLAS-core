#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_CONJ_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_CONJ_H 1

namespace ulmBLAS {

template <typename IndexType, typename VX>
    void
    conj(IndexType      n,
         VX             *x,
         IndexType      incX);

} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/conj.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_CONJ_H 1
