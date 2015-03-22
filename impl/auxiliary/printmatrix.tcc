/*
 * Copyright (C) 2014, The University of Texas at Austin
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of The University of Texas at Austin nor the names
 *    of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * Copyright (C) 2014-2015, Michael Lehn
 *
 * ulmBLAS adopted general ideas from BLIS.  Using micro kernels from BLIS
 * only requires minor modifications,
 *
 */

#ifndef ULMBLAS_IMPL_AUXILIARY_PRINTMATRIX_TCC
#define ULMBLAS_IMPL_AUXILIARY_PRINTMATRIX_TCC 1

#include <cstdio>
#include <complex>
#include <type_traits>
#include <ulmblas/impl/auxiliary/printmatrix.h>

namespace ulmBLAS {

template <typename T, typename IndexType>
void
printMatrix(IndexType m, IndexType n,
            const T *X, IndexType incRowX, IndexType incColX)
{
    if (std::is_same<double,T>::value) {
        for (IndexType i=0; i<m; ++i) {
            for (IndexType j=0; j<n; ++j) {
                //printf(" %7.4lf", X[i*incRowX+j*incColX]);
                printf(" %15.3lf", X[i*incRowX+j*incColX]);
            }
            printf("\n");
        }
        printf("\n");
    } else if (std::is_same<float, T>::value) {
        for (IndexType i=0; i<m; ++i) {
            for (IndexType j=0; j<n; ++j) {
                //printf(" %7.4lf", X[i*incRowX+j*incColX]);
                printf(" %15.3f", X[i*incRowX+j*incColX]);
            }
            printf("\n");
        }
        printf("\n");
    } else if (std::is_same<std::complex<double>, T>::value) {
        for (IndexType i=0; i<m; ++i) {
            for (IndexType j=0; j<n; ++j) {
                //printf(" %7.4lf", X[i*incRowX+j*incColX]);
                printf(" (%15.3lf, %15.3lf) ",
                       std::real(X[i*incRowX+j*incColX]),
                       std::imag(X[i*incRowX+j*incColX]));
            }
            printf("\n");
        }
        printf("\n");
    }
}

template <typename T, typename IndexType>
void
printSylMatrix(IndexType m,
               const T *X, IndexType incRowX, IndexType incColX)
{
    for (IndexType i=0; i<m; ++i) {
        for (IndexType j=0; j<m; ++j) {
            if (i>j) {
                printf(" %5.3lf", X[i*incRowX+j*incColX]);
            } else {
                printf(" %5.3lf", X[j*incRowX+i*incColX]);
            }
        }
        printf("\n");
    }
    printf("\n");
}

template <typename T, typename IndexType>
void
printSyuMatrix(IndexType m,
               const T *X, IndexType incRowX, IndexType incColX)
{
    for (IndexType i=0; i<m; ++i) {
        for (IndexType j=0; j<m; ++j) {
            if (i<j) {
                printf(" %5.3lf", X[i*incRowX+j*incColX]);
            } else {
                printf(" %5.3lf", X[j*incRowX+i*incColX]);
            }
        }
        printf("\n");
    }
    printf("\n");
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_AUXILIARY_PRINTMATRIX_TCC
