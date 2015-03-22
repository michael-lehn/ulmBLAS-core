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

#ifndef ULMBLAS_IMPL_LEVEL3_PACK_TRLPACK_TCC
#define ULMBLAS_IMPL_LEVEL3_PACK_TRLPACK_TCC 1

#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/level3/pack/trlpack.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>

namespace ulmBLAS {

template <typename IndexType, typename TL, typename Buffer>
static void
trlpack_MRxk(IndexType   k,
             bool        conj,
             bool        unit,
             const TL    *L,
             IndexType   incRowL,
             IndexType   incColL,
             Buffer      *buffer)
{
    const IndexType MR  = BlockSizeUGemm<Buffer>::MR;

    if (!conj) {
        for (IndexType j=0; j<k-MR; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                buffer[i] = L[i*incRowL];
            }
            buffer += MR;
            L      += incColL;
        }
    } else {
        for (IndexType j=0; j<k-MR; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                buffer[i] = conjugate(L[i*incRowL]);
            }
            buffer += MR;
            L      += incColL;
        }
    }
    for (IndexType j=0; j<MR; ++j) {
        for (IndexType i=0; i<j; ++i) {
            buffer[i] = Buffer(0);
        }
        buffer[j] = (unit) ? Buffer(1) : conjugate(L[j*incRowL], conj);
        for (IndexType i=j+1; i<MR; ++i) {
            buffer[i] = conjugate(L[i*incRowL], conj);
        }
        buffer += MR;
        L      += incColL;
    }
}

template <typename IndexType, typename TL, typename Buffer>
void
trlpack(IndexType   mc,
        bool        conj,
        bool        unit,
        const TL    *L,
        IndexType   incRowL,
        IndexType   incColL,
        Buffer      *buffer)
{
    const IndexType MR  = BlockSizeUGemm<Buffer>::MR;
    const IndexType mp  = mc / MR;
    const IndexType mr_ = mc % MR;

    for (IndexType i=0; i<mp; ++i) {
        trlpack_MRxk((i+1)*MR, conj, unit, L, incRowL, incColL, buffer);
        buffer += (i+1)*MR*MR;
        L      += MR*incRowL;
    }

    if (mr_>0) {
        for (IndexType j=0; j<mp*MR; ++j) {
            for (IndexType i=0; i<mr_; ++i) {
                buffer[i] = conjugate(L[i*incRowL], conj);
            }
            for (IndexType i=mr_; i<MR; ++i) {
                buffer[i] = Buffer(0);
            }
            buffer += MR;
            L      += incColL;
        }
        for (IndexType j=0; j<mr_; ++j) {
            for (IndexType i=0; i<j; ++i) {
                buffer[i] = Buffer(0);
            }
            buffer[j] = (unit) ? Buffer(1) : conjugate(L[j*incRowL], conj);
            for (IndexType i=j+1; i<mr_; ++i) {
                buffer[i] = conjugate(L[i*incRowL], conj);
            }
            for (IndexType i=mr_; i<MR; ++i) {
                buffer[i] = Buffer(0);
            }
            buffer += MR;
            L      += incColL;
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_PACK_TRLPACK_TCC
