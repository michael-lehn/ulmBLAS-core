#ifndef ULMBLAS_CXXBLAS_AUXILIARY_MACROS_H
#define ULMBLAS_CXXBLAS_AUXILIARY_MACROS_H 1

#ifndef CXXBLAS_DEBUG_OUT

#   ifdef CXXBLAS_DEBUG
#       include <iostream>
#       define CXXBLAS_DEBUG_OUT(x) std::cerr << "[ulmBLAS] " << x << std::endl;
#   else
#       define CXXBLAS_DEBUG_OUT(x)
#   endif

#endif

#endif // ULMBLAS_CXXBLAS_AUXILIARY_MACROS_H
