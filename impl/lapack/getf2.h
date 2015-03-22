#ifndef ULMBLAS_IMPL_LAPACK_GETF2_H
#define ULMBLAS_IMPL_LAPACK_GETF2_H 1

namespace ulmBLAS {

template <typename IndexType, typename T>
    IndexType
    getf2(IndexType    m,
          IndexType    n,
          T            *A,
          IndexType    incRowA,
          IndexType    incColA,
          IndexType    *piv,
          IndexType    incPiv);

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LAPACK_GETF2_H

#include <ulmblas/impl/lapack/getf2.tcc>
