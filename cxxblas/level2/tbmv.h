#ifndef ULMBLAS_CXXBLAS_LEVEL2_TBMV_H
#define ULMBLAS_CXXBLAS_LEVEL2_TBMV_H 1

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
    void
    tbmv(IndexType    n,
         IndexType    k,
         bool         colMajorA,
         bool         lowerA,
         bool         transA,
         bool         conjA,
         bool         unitDiagA,
         const TA     *A,
         IndexType    ldA,
         TX           *x,
         IndexType    incX);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/tbmv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_TBMV_H
