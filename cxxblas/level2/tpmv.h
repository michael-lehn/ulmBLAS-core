#ifndef ULMBLAS_CXXBLAS_LEVEL2_TPMV_H
#define ULMBLAS_CXXBLAS_LEVEL2_TPMV_H 1

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
    void
    tpmv(IndexType    n,
         bool         colMajorA,
         bool         lowerA,
         bool         transA,
         bool         conjA,
         bool         unitDiagA,
         const TA     *AP,
         TX           *x,
         IndexType    incX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/tpmv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_TPMV_H
