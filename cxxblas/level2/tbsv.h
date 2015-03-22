#ifndef ULMBLAS_CXXBLAS_LEVEL2_TBSV_H
#define ULMBLAS_CXXBLAS_LEVEL2_TBSV_H 1

namespace cxxblas {

template <typename IndexType, typename TA, typename TX>
    void
    tbsv(IndexType    n,
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

#include <ulmblas/cxxblas/level2/tbsv.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_TBSV_H
