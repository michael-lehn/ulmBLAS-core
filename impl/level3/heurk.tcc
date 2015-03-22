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

#ifndef ULMBLAS_IMPL_LEVEL3_HEURK_TCC
#define ULMBLAS_IMPL_LEVEL3_HEURK_TCC 1

#include <complex>
#include <ulmblas/impl/config/blocksize.h>
#include <ulmblas/impl/auxiliary/memorypool.h>
#include <ulmblas/impl/level1extensions/truscal.h>
#include <ulmblas/impl/level3/mkernel/mgemm.h>
#include <ulmblas/impl/level3/mkernel/msyurk.h>
#include <ulmblas/impl/level3/ukernel/ugemm.h>
#include <ulmblas/impl/level3/pack/gepack.h>
#include <ulmblas/impl/level3/heurk.h>

namespace ulmBLAS {

template <typename IndexType, typename Alpha, typename TA, typename Beta,
          typename TC>
void
heurk(IndexType    n,
      IndexType    k,
      const Alpha  &alpha,
      const TA     *A,
      IndexType    incRowA,
      IndexType    incColA,
      const Beta   &beta,
      TC           *C,
      IndexType    incRowC,
      IndexType    incColC)
{
    typedef decltype(Alpha(0)*TA(0))  T;

    const IndexType MC = BlockSize<T>::MC;

    const IndexType MR = BlockSizeUGemm<T>::MR;
    const IndexType NR = BlockSizeUGemm<T>::NR;

    const IndexType mb = (n+MC-1) / MC;
    const IndexType kb = (k+MC-1) / MC;

    const IndexType mc_ = n % MC;
    const IndexType kc_ = k % MC;

    static MemoryPool<T> memoryPool;

    if (n==0 || ((alpha==Alpha(0) || k==0) && beta==Beta(1))) {
        return;
    }

    if (alpha==Alpha(0) || k==0) {
        truscal(n, n, false, beta, C, incRowC, incColC);
        for (IndexType i=0; i<n; ++i) {
            C[i*(incRowC+incColC)] = std::real(C[i*(incRowC+incColC)]);
        }
        return;
    }

    T  *A_ = memoryPool.allocate(MC*MC+MR);
    T  *B_ = memoryPool.allocate(MC*MC+NR);

    for (IndexType j=0; j<mb; ++j) {
        IndexType nc = (j!=mb-1 || mc_==0) ? MC : mc_;

        for (IndexType l=0; l<kb; ++l) {
            IndexType kc    = (l!=kb-1 || kc_==0) ? MC   : kc_;
            Beta      beta_ = (l==0) ? beta : Beta(1);

            gepack_B(kc, nc, true,
                     &A[l*MC*incColA+j*MC*incRowA], incColA, incRowA,
                     B_);

            for (IndexType i=0; i<=j; ++i) {
                IndexType mc = (i!=mb-1 || mc_==0) ? MC : mc_;

                gepack_A(mc, kc, false,
                         &A[i*MC*incRowA+l*MC*incColA], incRowA, incColA,
                         A_);

                if (i==j) {
                    msyurk(mc, nc, kc, T(alpha), A_, B_, beta_,
                           &C[i*MC*incRowC+j*MC*incColC],
                           incRowC, incColC);
                } else {
                    mgemm(mc, nc, kc, T(alpha), A_, B_, beta_,
                          &C[i*MC*incRowC+j*MC*incColC],
                          incRowC, incColC);
                }
            }
        }
    }

    for (IndexType i=0; i<n; ++i) {
        C[i*(incRowC+incColC)] = std::real(C[i*(incRowC+incColC)]);
    }

    memoryPool.release(A_);
    memoryPool.release(B_);
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_HEURK_TCC
