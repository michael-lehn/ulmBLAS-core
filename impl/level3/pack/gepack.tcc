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
void
gepack_A(IndexType   mc,
         IndexType   kc,
         bool        conj,
         const TA    *A,
         IndexType   incRowA,
         IndexType   incColA,
         Buffer      *p)
{
    IndexType MR = BlockSize<Buffer>::MR;
    IndexType mp = (mc+MR-1) / MR;

    for (IndexType j=0; j<kc; ++j) {
        for (IndexType l=0; l<mp; ++l) {
            for (IndexType i0=0; i0<MR; ++i0) {
                IndexType i  = l*MR + i0;
                IndexType nu = l*MR*kc + j*MR + i0;
                p[nu]        = (i<mc) ? conjugate(A[i*incRowA+j*incColA],conj)
                                      : Buffer(0);
            }
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
         Buffer      *p)
{
    const IndexType NR = BlockSize<Buffer>::NR;
    const IndexType np = (nc+NR-1) / NR;

    for (IndexType l=0; l<np; ++l) {
        for (IndexType i=0; i<kc; ++i) {
            for (IndexType j0=0; j0<NR; ++j0) {
                IndexType j  = l*NR+j0;
                IndexType nu = l*NR*kc + i*NR + j0;
                p[nu]        = (j<nc) ? conjugate(B[i*incRowB+j*incColB],conj)
                                      : Buffer(0);
            }
        }
    }
}

template <typename IndexType, typename Buffer, typename TB>
void
geunpack_B(IndexType    kc,
           IndexType    nc,
           bool         conj,
           const Buffer *p,
           TB           *B,
           IndexType    incRowB,
           IndexType    incColB)
{
    IndexType NR = BlockSize<Buffer>::NR;
    IndexType np = (nc+NR-1) / NR;

    for (IndexType l=0; l<np; ++l) {
        for (IndexType i=0; i<kc; ++i) {
            for (IndexType j0=0; j0<NR; ++j0) {
                IndexType j  = l*NR+j0;
                IndexType nu = l*NR*kc + i*NR + j0;
                if (j<nc) {
                    B[i*incRowB+j*incColB] = p[nu];
                }
            }
        }
    }
}


} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_PACK_GEPACK_TCC
