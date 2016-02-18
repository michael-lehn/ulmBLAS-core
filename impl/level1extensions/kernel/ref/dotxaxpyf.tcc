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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOTXAXPYF_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOTXAXPYF_TCC 1

#include <type_traits>
#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/config/fusefactor.h>
#include <ulmblas/impl/level1extensions/kernel/ref/dotxaxpyf.h>

#ifdef DDOTXAXPYF_FUSEFACTOR
#undef DDOTXAXPYF_FUSEFACTOR
#endif

#define DDOTXAXPYF_FUSEFACTOR  2

namespace ulmBLAS { namespace ref {

template <typename IndexType, typename Alpha, typename VA, typename MX,
          typename VY, typename VZ, typename Rho>
void
dotxaxpyf(IndexType      n,
          bool           conjX,
          bool           conjXt,
          bool           conjY,
          const Alpha    &alpha,
          const VA       *a,
          IndexType      incA,
          const MX       *X,
          IndexType      incRowX,
          IndexType      incColX,
          const VY       *y,
          IndexType      incY,
          VZ             *z,
          IndexType      incZ,
          Rho            *rho,
          IndexType      incRho)
{
    typedef decltype(Alpha(0)*VA(0)*MX(0)*VY(0)*VZ(0))  T;

    const IndexType bf = FuseFactor<T>::dotxaxpyf;

    for (IndexType l=0; l<bf; ++l) {
        rho[l*incRho] = 0;
    }

    if (n<=0) {
        return;
    }

    for (IndexType i=0; i<n; ++i) {
        const VY y_i = conjugate(y[i*incY], conjY);

        for (IndexType l=0; l<bf; ++l) {
            const MX x_il  = conjugate(X[i*incRowX+l*incColX], conjX);
            const MX xt_il = conjugate(X[i*incRowX+l*incColX], conjXt);

            rho[l*incRho] += xt_il * y_i;
            z[i*incZ]     += alpha*a[l*incA]*x_il;
        }
    }
}

} } // namespace ref, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOTXAXPYF_TCC 1
