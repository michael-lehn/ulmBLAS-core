#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_CONJ_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_CONJ_TCC 1

#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/level1extensions/conj.h>

namespace ulmBLAS {

template <typename IndexType, typename VX>
void
conj(IndexType      n,
     VX             *x,
     IndexType      incX)
{
    for (IndexType i=0; i<n; ++i) {
        x[i*incX] = conjugate(x[i*incX]);
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_CONJ_TCC 1
