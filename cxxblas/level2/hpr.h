#ifndef ULMBLAS_CXXBLAS_LEVEL2_HPR_H
#define ULMBLAS_CXXBLAS_LEVEL2_HPR_H 1

namespace cxxblas {

template <typename IndexType, typename Alpha, typename TX, typename TA>
    void
    hpr(IndexType    n,
        const Alpha  &alpha,
        const TX     *x,
        IndexType    incX,
        bool         colMajorA,
        bool         lowerA,
        TA           *A);

} // namespace cxxblas

#include <ulmblas/cxxblas/level2/hpr.tcc>

#endif // ULMBLAS_CXXBLAS_LEVEL2_HPR_H
