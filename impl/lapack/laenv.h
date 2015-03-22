#ifndef ULMBLAS_IMPL_LAPACK_LAENV_H
#define ULMBLAS_IMPL_LAPACK_LAENV_H 1

namespace ulmBLAS {

template <typename T>
    long
    laenv(long        spec,
          const char  *name,
          const char  *opts,
          long        n1 = -1,
          long        n2 = -1,
          long        n3 = -1,
          long        n4 = -1);

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LAPACK_LAENV_H

#include <ulmblas/impl/lapack/laenv.tcc>
