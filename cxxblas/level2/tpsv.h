#ifndef ULMBLAS_CXXBLAS_LEVEL2_TPSV_H
#define ULMBLAS_CXXBLAS_LEVEL2_TPSV_H 1

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
    void
    tpsv(IndexType    n,
         bool         colMajorA,
         bool         lowerA,
         bool         transA,
         bool         conjA,
         bool         unitDiagA,
         const TA     *AP,
         TX           *x,
         IndexType    incX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/tpsv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_TPSV_H
