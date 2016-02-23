#ifndef ULMBLAS_EXTERNAL_BLIS_SCAL_H
#define ULMBLAS_EXTERNAL_BLIS_SCAL_H 1

namespace cxxblas {

template <typename IndexType, typename T>
    void
    scal(IndexType      n,
         const T        &alpha,
         T              *x,
         IndexType      incX);

} // namespace cxxblas

#endif // ULMBLAS_EXTERNAL_BLIS_SCAL_H

#include <ulmblas/external/blis/scal.tcc>
