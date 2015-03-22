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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_AXPYF_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_AXPYF_TCC 1

#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/level1extensions/kernel/ref/axpyf.h>

namespace ulmBLAS { namespace ref {

template <typename IndexType, typename Alpha, typename VA, typename VX,
           typename VY>
typename std::enable_if<std::is_integral<IndexType>::value
                  && FuseFactor<decltype(Alpha(0)*VA(0)*VX(0)+VY(0))>::axpyf==4,
    void>::type
axpyf(IndexType      n,
      const Alpha    &alpha,
      const VA       *a,
      IndexType      incA,
      const VX       *X,
      IndexType      incRowX,
      IndexType      incColX,
      VY             *y,
      IndexType      incY)
{
    const VX *x0 = &X[0*incColX];
    const VX *x1 = &X[1*incColX];
    const VX *x2 = &X[2*incColX];
    const VX *x3 = &X[3*incColX];

    const VA &a0 = a[0*incA];
    const VA &a1 = a[1*incA];
    const VA &a2 = a[2*incA];
    const VA &a3 = a[3*incA];

    for (IndexType i=0; i<n; ++i) {
        y[i*incY] += alpha*(a0*x0[i*incRowX]+a1*x1[i*incRowX]
                           +a2*x2[i*incRowX]+a3*x3[i*incRowX]);
    }
}

template <typename IndexType, typename Alpha, typename VA, typename VX,
           typename VY>
typename std::enable_if<std::is_integral<IndexType>::value
                  && FuseFactor<decltype(Alpha(0)*VA(0)*VX(0)+VY(0))>::axpyf!=4,
    void>::type
axpyf(IndexType      n,
      const Alpha    &alpha,
      const VA       *a,
      IndexType      incA,
      const VX       *X,
      IndexType      incRowX,
      IndexType      incColX,
      VY             *y,
      IndexType      incY)
{
    typedef decltype(Alpha(0)*VA(0)*VX(0)+VY(0))    T;

    const IndexType ff = FuseFactor<T>::axpyf;

    for (IndexType i=0; i<n; ++i) {
        for (IndexType j=0; j<ff; ++j) {
            y[i*incY] += alpha*a[j*incA]*X[i*incRowX+j*incColX];
        }
    }
}


template <typename IndexType, typename Alpha, typename VA, typename VX,
           typename VY>
typename std::enable_if<std::is_integral<IndexType>::value
                 && FuseFactor<decltype(Alpha(0)*VA(0)*VX(0)+VY(0))>::acxpyf==4,
    void>::type
acxpyf(IndexType      n,
       const Alpha    &alpha,
       const VA       *a,
       IndexType      incA,
       const VX       *X,
       IndexType      incRowX,
       IndexType      incColX,
       VY             *y,
       IndexType      incY)
{
    const VX *x0 = &X[0*incColX];
    const VX *x1 = &X[1*incColX];
    const VX *x2 = &X[2*incColX];
    const VX *x3 = &X[3*incColX];

    const VA &a0 = a[0*incA];
    const VA &a1 = a[1*incA];
    const VA &a2 = a[2*incA];
    const VA &a3 = a[3*incA];

    for (IndexType i=0; i<n; ++i) {
        y[i*incY] += alpha*(a0*conjugate(x0[i*incRowX])
                           +a1*conjugate(x1[i*incRowX])
                           +a2*conjugate(x2[i*incRowX])
                           +a3*conjugate(x3[i*incRowX]));
    }
}

template <typename IndexType, typename Alpha, typename VA, typename VX,
           typename VY>
typename std::enable_if<std::is_integral<IndexType>::value
                 && FuseFactor<decltype(Alpha(0)*VA(0)*VX(0)+VY(0))>::acxpyf!=4,
    void>::type
acxpyf(IndexType      n,
       const Alpha    &alpha,
       const VA       *a,
       IndexType      incA,
       const VX       *X,
       IndexType      incRowX,
       IndexType      incColX,
       VY             *y,
       IndexType      incY)
{
    typedef decltype(Alpha(0)*VA(0)*VX(0)+VY(0))    T;

    const IndexType ff = FuseFactor<T>::axpyf;

    for (IndexType i=0; i<n; ++i) {
        for (IndexType j=0; j<ff; ++j) {
            y[i*incY] += alpha*a[j*incA]*conjugate(X[i*incRowX+j*incColX]);
        }
    }
}

} } // namespace ref, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_AXPYF_TCC
