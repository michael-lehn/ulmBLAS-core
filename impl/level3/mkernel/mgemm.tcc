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

#ifndef ULMBLAS_IMPL_LEVEL3_MKERNEL_MGEMM_TCC
#define ULMBLAS_IMPL_LEVEL3_MKERNEL_MGEMM_TCC 1

#include <ulmblas/impl/level3/mkernel/mgemm.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>

namespace ulmBLAS {


template <typename IndexType, typename T, typename Beta, typename TC>
void
mgemm(IndexType     mc,
      IndexType     nc,
      IndexType     kc,
      const T       &alpha,
      const T       *A_,
      const T       *B_,
      const Beta    &beta,
      TC            *C,
      IndexType     incRowC,
      IndexType     incColC)
{
    const IndexType MR = BlockSize<T>::MR;
    const IndexType NR = BlockSize<T>::NR;

    const IndexType mp = (mc+MR-1) / MR;
    const IndexType np = (nc+NR-1) / NR;

    const IndexType mr_ = mc % MR;
    const IndexType nr_ = nc % NR;

    IndexType mr, nr;

    const T *nextA;
    const T *nextB;

    for (IndexType j=0; j<np; ++j) {
        nr    = (j!=np-1 || nr_==0) ? NR : nr_;
        nextB = &B_[j*kc*NR];

        for (IndexType i=0; i<mp; ++i) {
            mr    = (i!=mp-1 || mr_==0) ? MR : mr_;
            nextA = &A_[(i+1)*kc*MR];

            if (i==mp-1) {
                nextA = A_;
                nextB = &B_[(j+1)*kc*NR];
                if (j==np-1) {
                    nextB = B_;
                }
            }

            if (mr==MR && nr==NR) {
                ugemm(kc,
                      alpha, &A_[i*kc*MR], &B_[j*kc*NR],
                      beta,
                      &C[i*MR*incRowC+j*NR*incColC], incRowC, incColC,
                      nextA, nextB);
            } else {
                // Call the buffered micro kernel
                ugemm(mr, nr, kc,
                      alpha, &A_[i*kc*MR], &B_[j*kc*NR],
                      beta,
                      &C[i*MR*incRowC+j*NR*incColC], incRowC, incColC,
                      nextA, nextB);
            }
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_MKERNEL_MGEMM_TCC
