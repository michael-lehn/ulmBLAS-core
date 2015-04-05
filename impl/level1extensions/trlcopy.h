#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLCOPY_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLCOPY_H 1

namespace ulmBLAS {

template <typename IndexType, typename MX, typename MY>
    void
    trlcopy(IndexType    m,
            IndexType    n,
            bool         unit,
            bool         conjA,
            MX           *X,
            IndexType    incRowX,
            IndexType    incColX,
            MY           *Y,
            IndexType    incRowY,
            IndexType    incColY);

} // namespace ulmBLAS

#include <ulmblas/impl/level1extensions/trlcopy.tcc>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLCOPY_H 1
