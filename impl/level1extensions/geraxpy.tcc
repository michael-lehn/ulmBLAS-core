/*
 * Copyright (C) 2014, The University of Texas at Austin
 * Copyright (C) 2014-2015, Michael Lehn
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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GERAXPY_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GERAXPY_TCC 1

#include <ulmblas/impl/level1extensions/geaxpy.h>
#include <ulmblas/impl/level1/axpy.h>
#include <ulmblas/impl/auxiliary/conjugate.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename MX, typename MY>
void
geraxpy(IndexType      m,
        IndexType      n,
        const Alpha    &alpha_,
        const MX       *X,
        IndexType      incRowX,
        IndexType      incColX,
        MY             *Y,
        IndexType      incRowY,
        IndexType      incColY)
{
    assert(alpha_!=Alphe(0));
    const IndexType    UnitStride(1);
    const Alpha        alpha = Alpha(1)/alpha_;

    if (m<=0 || n<=0 || alpha==Alpha(0)) {
        return;
    }

    if (incRowX==UnitStride && incRowY==UnitStride) {
//
//      X and Y are both column major
//
        for (IndexType j=0; j<n; ++j) {
            axpy(m, alpha,
                 &X[j*incColX], UnitStride,
                 &Y[j*incColY], UnitStride);
        }
    } else if (incColX==UnitStride && incColY==UnitStride) {
//
//      X and Y are both row major
//
        for (IndexType i=0; i<m; ++i) {
            axpy(n, alpha,
                 &X[i*incRowX], UnitStride,
                 &Y[i*incRowY], UnitStride);
        }
    } else {
//
//      General case
//
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<m; ++i) {
                Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
            }
        }
    }
}

template <typename IndexType, typename Alpha, typename MX, typename MY>
void
geracxpy(IndexType      m,
         IndexType      n,
         const Alpha    &alpha_,
         const MX       *X,
         IndexType      incRowX,
         IndexType      incColX,
         MY             *Y,
         IndexType      incRowY,
         IndexType      incColY)
{
    assert(alpha_!=Alphe(0));
    const IndexType    UnitStride(1);
    const Alpha        alpha = Alpha(1)/alpha_;

    if (m<=0 || n<=0 || alpha==Alpha(0)) {
        return;
    }

    if (incRowX==UnitStride && incRowY==UnitStride) {
//
//      X and Y are both column major
//
        for (IndexType j=0; j<n; ++j) {
            acxpy(m, alpha,
                  &X[j*incColX], UnitStride,
                  &Y[j*incColY], UnitStride);
        }
    } else if (incColX==UnitStride && incColY==UnitStride) {
//
//      X and Y are both row major
//
        for (IndexType i=0; i<m; ++i) {
            acxpy(n, alpha,
                  &X[i*incRowX], UnitStride,
                  &Y[i*incRowY], UnitStride);
        }
    } else {
//
//      General case
//
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<m; ++i) {
                Y[i*incRowY+j*incColY] += alpha
                                          *conjugate(X[i*incRowX+j*incColX]);
            }
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GERAXPY_TCC 1
