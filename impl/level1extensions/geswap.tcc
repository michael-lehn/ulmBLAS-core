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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_TCC 1

#include <ulmblas/impl/level1extensions/geswap.h>

namespace ulmBLAS {

template <typename IndexType, typename MX, typename MY>
void
geswap(IndexType      m,
       IndexType      n,
       const MX       *X,
       IndexType      incRowX,
       IndexType      incColX,
       MY             *Y,
       IndexType      incRowY,
       IndexType      incColY)
{
    if (incRowX<incColX) {
        for (IndexType j=0; j<n; ++j) {
            swap(m, &X[j*incColX], incRowX, &Y[j*incColY], incRowY);
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            swap(n, &X[i*incRowX], incColX, &Y[i*incRowY], incColY);
        }
    }
}

template <typename IndexType, typename MX, typename MY>
void
geswapc(IndexType      m,
        IndexType      n,
        const MX       *X,
        IndexType      incRowX,
        IndexType      incColX,
        MY             *Y,
        IndexType      incRowY,
        IndexType      incColY)
{
    if (incRowX<incColX) {
        for (IndexType j=0; j<n; ++j) {
            swapc(m, &X[j*incColX], incRowX, &Y[j*incColY], incRowY);
        }
    } else {
        for (IndexType i=0; i<m; ++i) {
            swapc(n, &X[i*incRowX], incColX, &Y[i*incRowY], incColY);
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_GESWAP_TCC 1
