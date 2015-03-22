#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESCAL_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESCAL_H 1

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename MA>
    void
    gescal(IndexType    m,
           IndexType    n,
           const Alpha  &alpha,
           MA           *A,
           IndexType    incRowA,
           IndexType    incColA);

} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/gescal.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESCAL_H 1
