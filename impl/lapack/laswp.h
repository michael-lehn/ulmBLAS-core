#ifndef ULMBLAS_IMPL_LAPACK_LASWP_H
#define ULMBLAS_IMPL_LAPACK_LASWP_H 1

namespace ulmBLAS {

template <typename IndexType, typename T>
    void
    laswp(IndexType    n,
          T            *A,
          IndexType    incRowA,
          IndexType    incColA,
          IndexType    k1,
          IndexType    k2,
          IndexType    *piv,
          IndexType    incPiv);

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LAPACK_LASWP_H

#include <ulmblas/impl/lapack/laswp.tcc>
