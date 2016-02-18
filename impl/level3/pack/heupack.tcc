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

#ifndef ULMBLAS_IMPL_LEVEL3_PACK_HEUPACK_TCC
#define ULMBLAS_IMPL_LEVEL3_PACK_HEUPACK_TCC 1

#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/level3/pack/helpack.h>
#include <ulmblas/impl/level3/pack/gepack.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>

namespace ulmBLAS {

template <typename IndexType, typename TA, typename Buffer>
static void
heupack_mrxmr(IndexType   mr,
              const TA    *A,
              IndexType   incRowA,
              IndexType   incColA,
              Buffer      *buffer)
{
    const IndexType MR  = BlockSizeUGemm<Buffer>::MR;

    for (IndexType j=0; j<mr; ++j) {
        for (IndexType i=0; i<mr; ++i) {
            buffer[i] = (i<j)  ? A[i*incRowA+j*incColA]
                      : (i==j) ? std::real(A[i*incRowA+j*incColA])
                      : conjugate(A[j*incRowA+i*incColA]);
        }
        for (IndexType i=mr; i<MR; ++i) {
            buffer[i] = Buffer(0);
        }
        buffer += MR;
    }
}

template <typename IndexType, typename TA, typename Buffer>
void
heupack(IndexType   mc,
        const TA    *A,
        IndexType   incRowA,
        IndexType   incColA,
        Buffer      *buffer)
{
    const IndexType MR  = BlockSizeUGemm<Buffer>::MR;
    const IndexType mp  = mc / MR;
    const IndexType mr_ = mc % MR;

    for (IndexType i=0; i<mp; ++i) {
        gepack_A(MR, i*MR, true,
                 &A[i*MR*incColA], incColA, incRowA,
                 buffer);
        buffer += MR*i*MR;

        heupack_mrxmr(MR,
                      &A[i*MR*incRowA+i*MR*incColA], incRowA, incColA,
                      buffer);
        buffer += MR*MR;

        gepack_A(MR, mc-(i+1)*MR, false,
                 &A[i*MR*incRowA+(i+1)*MR*incColA], incRowA, incColA,
                 buffer);
        buffer += MR*(mc-(i+1)*MR);
    }

    if (mr_>0) {
        gepack_A(mr_, mc-mr_, true,
                 &A[mp*MR*incColA], incColA, incRowA,
                 buffer);
        buffer += MR*(mc-mr_);

        heupack_mrxmr(mr_,
                      &A[mp*MR*incRowA+mp*MR*incColA], incRowA, incColA,
                      buffer);
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_PACK_HEUPACK_TCC
