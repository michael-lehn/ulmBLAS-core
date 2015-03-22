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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_REF_UGEMM_TCC
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_REF_UGEMM_TCC 1

#include <ulmblas/impl/level3/ukernel/ref/ugemm.h>

namespace ulmBLAS { namespace ref {

template <typename IndexType, typename T>
void
ugemm(IndexType   kc,
      const T     &alpha,
      const T     *A,
      const T     *B,
      const T     &beta,
      T           *C,
      IndexType   incRowC,
      IndexType   incColC,
      const T     *,
      const T     *)
{
    const IndexType MR = BlockSizeUGemm<T>::MR;
    const IndexType NR = BlockSizeUGemm<T>::NR;

    T AB[MR*NR];

    for (IndexType i=0; i<MR*NR; ++i) {
        AB[i] = T(0);
    }

    for (IndexType l=0; l<kc; ++l) {
        for (IndexType j=0; j<NR; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                AB[i+j*MR] += A[i]*B[j];
            }
        }
        A += MR;
        B += NR;
    }

    if (beta==T(0)) {
        for (IndexType j=0; j<NR; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                C[i*incRowC+j*incColC] = T(0);
            }
        }
    } else {
        for (IndexType j=0; j<NR; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                C[i*incRowC+j*incColC] *= beta;
            }
        }
    }

    if (alpha==T(1)) {
        for (IndexType j=0; j<NR; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                C[i*incRowC+j*incColC] += AB[i+j*MR];
            }
        }
    } else {
        for (IndexType j=0; j<NR; ++j) {
            for (IndexType i=0; i<MR; ++i) {
                C[i*incRowC+j*incColC] += alpha*AB[i+j*MR];
            }
        }
    }
}

} } // namespace ref, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_REF_UGEMM_TCC 1
