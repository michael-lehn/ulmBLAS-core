#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLSCAL_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLSCAL_H 1

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename MA>
    void
    trlscal(IndexType    m,
            IndexType    n,
            bool         unit,
            const Alpha  &alpha,
            MA           *A,
            IndexType    incRowA,
            IndexType    incColA);

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLSCAL_H 1

#include <ulmblas/impl/level1extensions/trlscal.tcc>
