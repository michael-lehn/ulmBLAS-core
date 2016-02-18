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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_TCC 1

#include <ulmblas/impl/level1extensions/gecotr.h>
#include <ulmblas/impl/level1/swap.h>
#include <ulmblas/impl/auxiliary/conjugate.h>

namespace ulmBLAS {

template <typename IndexType, typename MX>
void
gecotr(IndexType      n,
       bool           transX,
       bool           conjX,
       MX             *X,
       IndexType      incRowX,
       IndexType      incColX)
{
    if (!transX && !conjX) {
        return;
    }

    if (transX) {
        if (!conjX) {
            for (IndexType i=0; i<n; ++i) {
                swap(n-1-i,
                     &X[(i+1)*incRowX+i*incColX], incRowX,
                     &X[i*incRowX+(i+1)*incColX], incColX);
            }
        } else {
            for (IndexType i=0; i<n; ++i) {
                swapc(n-1-i,
                      &X[(i+1)*incRowX+i*incColX], incRowX,
                      &X[i*incRowX+(i+1)*incColX], incColX);
                X[i*(incRowX+incColX)] = conjugate(X[i*(incRowX+incColX)]);
            }
        }
    } else {
        for (IndexType j=0; j<n; ++j) {
            for (IndexType i=0; i<n; ++i) {
                X[i*incRowX+j*incColX] = conjugate(X[i*incRowX+j*incColX]);
            }
        }
    }
}


} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GECOTR_TCC 1
