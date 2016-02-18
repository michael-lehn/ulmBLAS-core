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

#ifndef ULMBLAS_IMPL_LEVEL3_PACK_GEPACK_TCC
#define ULMBLAS_IMPL_LEVEL3_PACK_GEPACK_TCC 1

#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/level3/pack/gepack.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>

namespace ulmBLAS {

template <typename IndexType, typename TA, typename Buffer>
static void
pack_MRxk(IndexType   k,
          bool        conj,
          const TA    *A,
          IndexType   incRowA,
          IndexType   incColA,
          Buffer      *buffer)
{
    const IndexType MR = BlockSize<Buffer>::MR;

    if (!conj) {
        for (IndexType j=0; j<k; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                buffer[i] = A[i*incRowA];
            }
            buffer += MR;
            A      += incColA;
        }
    } else {
        for (IndexType j=0; j<k; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                buffer[i] = conjugate(A[i*incRowA]);
            }
            buffer += MR;
            A      += incColA;
        }
    }
}

template <typename IndexType, typename TA, typename Buffer>
void
gepack_A(IndexType   mc,
         IndexType   kc,
         bool        conj,
         const TA    *A,
         IndexType   incRowA,
         IndexType   incColA,
         Buffer      *buffer)
{
    const IndexType MR  = BlockSize<Buffer>::MR;
    const IndexType mp  = mc / MR;
    const IndexType mr_ = mc % MR;

    for (IndexType i=0; i<mp; ++i) {
        pack_MRxk(kc, conj, A, incRowA, incColA, buffer);
        buffer += kc*MR;
        A      += MR*incRowA;
    }
    if (mr_>0) {
        for (IndexType j=0; j<kc; ++j) {
            for (IndexType i=0; i<mr_; ++i) {
                buffer[i] = (!conj) ? A[i*incRowA] : conjugate(A[i*incRowA]);
            }
            for (IndexType i=mr_; i<MR; ++i) {
                buffer[i] = Buffer(0);
            }
            buffer += MR;
            A      += incColA;
        }
    }
}

template <typename IndexType, typename TB, typename Buffer>
static void
pack_kxNR(IndexType   k,
          bool        conj,
          const TB    *B,
          IndexType   incRowB,
          IndexType   incColB,
          Buffer      *buffer)
{
    const IndexType NR = BlockSize<Buffer>::NR;

    if (!conj) {
        for (IndexType i=0; i<k; ++i) {
            for (IndexType j=0; j<NR; ++j) {
                buffer[j] = B[j*incColB];
            }
            buffer += NR;
            B      += incRowB;
        }
    } else {
        for (IndexType i=0; i<k; ++i) {
            for (IndexType j=0; j<NR; ++j) {
                buffer[j] = conjugate(B[j*incColB]);
            }
            buffer += NR;
            B      += incRowB;
        }
    }
}

template <typename IndexType, typename TB, typename Buffer>
void
gepack_B(IndexType   kc,
         IndexType   nc,
         bool        conj,
         const TB    *B,
         IndexType   incRowB,
         IndexType   incColB,
         Buffer      *buffer)
{
    const IndexType NR  = BlockSize<Buffer>::NR;
    const IndexType np  = nc / NR;
    const IndexType nr_ = nc % NR;

    for (IndexType j=0; j<np; ++j) {
        pack_kxNR(kc, conj, B, incRowB, incColB, buffer);
        buffer += kc*NR;
        B      += NR*incColB;
    }
    if (nr_>0) {
        for (IndexType i=0; i<kc; ++i) {
            for (IndexType j=0; j<nr_; ++j) {
                buffer[j] = (!conj) ? B[j*incColB] : conjugate(B[j*incColB]);
            }
            for (IndexType j=nr_; j<NR; ++j) {
                buffer[j] = Buffer(0);
            }
            buffer += NR;
            B      += incRowB;
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_PACK_GEPACK_TCC
