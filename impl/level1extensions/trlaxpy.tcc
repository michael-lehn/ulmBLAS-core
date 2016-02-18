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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLAXPY_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLAXPY_TCC 1

#include <algorithm>
#include <ulmblas/impl/level1/axpy.h>
#include <ulmblas/impl/level1extensions/trlaxpy.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename MX, typename MY>
void
trlaxpy(IndexType    m,
        IndexType    n,
        bool         unit,
        const Alpha  &alpha,
        MX           *X,
        IndexType    incRowX,
        IndexType    incColX,
        MY           *Y,
        IndexType    incRowY,
        IndexType    incColY)
{
    const IndexType    UnitStride(1);

    if (m<=0 || n<=0 || alpha==Alpha(0)) {
        return;
    }

    if (unit) {
        trlaxpy(m-1, n-1, false, alpha,
                &X[1*incRowX], incRowX, incColX,
                &Y[1*incRowY], incRowY, incColY);
        return;
    }

    if (incRowX==UnitStride && incRowY==UnitStride) {
        const IndexType k = std::min(m, n);
        for (IndexType j=0; j<k; ++j) {
            axpy(m-j, alpha,
                 &X[j*(incRowX+incColX)], UnitStride,
                 &Y[j*(incRowY+incColY)], UnitStride);
        }
    } else if (incColX==UnitStride && incColY==UnitStride) {
        for (IndexType i=0; i<m; ++i) {
            axpy(std::min(i+1,n), alpha,
                 &X[i*incRowX], UnitStride,
                 &Y[i*incRowY], UnitStride);
        }
    } else {
        const IndexType k = std::min(m, n);
        for (IndexType j=0; j<k; ++j) {
            for (IndexType i=j; i<m; ++i) {
                Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
            }
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_TRLAXPY_TCC 1
