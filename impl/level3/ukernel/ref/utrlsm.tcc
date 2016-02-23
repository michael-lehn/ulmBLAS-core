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

#ifndef ULMBLAS_IMPL_LEVEL3_UKERNEL_REF_UTRLSM_TCC
#define ULMBLAS_IMPL_LEVEL3_UKERNEL_REF_UTRLSM_TCC 1

#include <complex>
#include <ulmblas/impl/level3/ukernel/ugemm.h>
#include <ulmblas/impl/level1extensions/gecopy.h>
#include <ulmblas/impl/level3/ukernel/ref/utrlsm.h>

namespace ulmBLAS { namespace ref {

template <typename T>
typename std::enable_if<BlockSize<T>::vlen == 0,
         void>::type
utrlsm(const T *A, T *B)
{
    typedef std::size_t IndexType;

    const IndexType MR = BlockSize<T>::MR;
    const IndexType NR = BlockSize<T>::NR;

    T   C_[MR*NR];

    for (IndexType i=0; i<MR; ++i) {
        for (IndexType j=0; j<NR; ++j) {
            C_[i+j*MR] = B[i*NR+j];
        }
    }

    for (IndexType i=0; i<MR; ++i) {
        for (IndexType j=0; j<NR; ++j) {
            C_[i+j*MR] *= A[i];
            for (IndexType l=i+1; l<MR; ++l) {
                C_[l+j*MR] -= A[l]*C_[i+j*MR];
            }
        }
        A += MR;
    }
    for (IndexType i=0; i<MR; ++i) {
        for (IndexType j=0; j<NR; ++j) {
            B[i*NR+j] = C_[i+j*MR];
        }
   }
}

template <typename T>
typename std::enable_if<BlockSize<T>::vlen != 0,
         void>::type
utrlsm(const T *A, T *B)
{
    typedef std::size_t IndexType;
    typedef T           vx __attribute__((vector_size(BlockSize<T>::rwidth/8)));

    static constexpr IndexType vlen = BlockSize<T>::vlen;
    static constexpr IndexType MR   = BlockSize<T>::MR;
    static constexpr IndexType NR   = BlockSize<T>::NR/vlen;

    A = (const T*) __builtin_assume_aligned (A, BlockSize<T>::align);
    B = ( T*)      __builtin_assume_aligned (B, BlockSize<T>::align);

    vx C_[MR*NR];
    vx *B_ = (vx *)B;

    for (IndexType i=0; i<MR*NR; ++i) {
        C_[i] = B_[i];
    }

    for (IndexType i=0; i<MR; ++i) {
        for (IndexType j=0; j<NR; ++j) {
            C_[i*NR+j] *= A[i];
        }
        for (IndexType l=i+1; l<MR; ++l) {
            for (IndexType j=0; j<NR; ++j) {
                C_[l*NR+j] -= A[l]*C_[i*NR+j];
            }
        }
        A += MR;
    }
    for (IndexType i=0; i<MR*NR; ++i) {
        B_[i] = C_[i];
    }
}


template <typename T>
void
utrlsm(const std::complex<T> *A, std::complex<T> *B)
{
    typedef std::size_t IndexType;

    const IndexType MR = BlockSize<std::complex<T>>::MR;
    const IndexType NR = BlockSize<std::complex<T>>::NR;

    T   C_[2*MR*NR];

    const T *rA  = reinterpret_cast<const T*>(A);
    const T *iA  = rA+1;

    T       *rB  = reinterpret_cast<T*>(B);
    T       *iB  = rB+1;

    T       *rC_ = reinterpret_cast<T*>(C_);
    T       *iC_ = rC_+1;

    for (IndexType i=0; i<MR; ++i) {
        for (IndexType j=0; j<NR; ++j) {
            rC_[(i+j*MR)*2] = rB[(i*NR+j)*2];
            iC_[(i+j*MR)*2] = iB[(i*NR+j)*2];
        }
    }

    for (IndexType i=0; i<MR; ++i) {
        for (IndexType j=0; j<NR; ++j) {
            T real          = rC_[2*(i+j*MR)]*rA[2*i] - iC_[2*(i+j*MR)]*iA[2*i];
            iC_[(i+j*MR)*2] = rC_[2*(i+j*MR)]*iA[2*i] + iC_[2*(i+j*MR)]*rA[2*i];
            rC_[(i+j*MR)*2] = real;

            for (IndexType l=i+1; l<MR; ++l) {
                rC_[2*(l+j*MR)] -= rC_[2*(i+j*MR)]*rA[2*l] - iC_[2*(i+j*MR)]*iA[2*l];
                iC_[2*(l+j*MR)] -= rC_[2*(i+j*MR)]*iA[2*l] + iC_[2*(i+j*MR)]*rA[2*l];
            }
        }
        rA += 2*MR;
        iA += 2*MR;
    }
    for (IndexType i=0; i<MR; ++i) {
        for (IndexType j=0; j<NR; ++j) {
            rB[(i*NR+j)*2] = rC_[(i+j*MR)*2];
            iB[(i*NR+j)*2] = iC_[(i+j*MR)*2];
        }
   }
}

} } // namespace ref, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL3_UKERNEL_REF_UTRLSM_TCC 1
