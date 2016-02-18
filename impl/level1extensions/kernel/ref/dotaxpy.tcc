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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOTAXPY_TCC
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOTAXPY_TCC 1

#include <complex>
#include <ulmblas/impl/auxiliary/conjugate.h>
#include <ulmblas/impl/level1extensions/kernel/ref/dotaxpy.h>

namespace ulmBLAS { namespace ref {

//
//  Fuse the computations:
//  (1) $z \leftarrow z + \alpha x$
//  (2) $\rho = x^T y$
//
//  Arguments $x$, $x^T$ or $y$ can be conjugated in the computation
//
template <typename IndexType, typename Alpha, typename VX, typename VY,
          typename VZ, typename Rho>
void
dotaxpy(IndexType      n,
        bool           conjX,
        bool           conjXt,
        bool           conjY,
        const Alpha    &alpha,
        const VX       *x,
        IndexType      incX,
        const VY       *y,
        IndexType      incY,
        VZ             *z,
        IndexType      incZ,
        Rho            &rho)
{
    rho = Rho(0);
    for (IndexType i=0; i<n; ++i) {
        z[i*incZ] += alpha*conjugate(x[i*incX], conjX);
        rho       += conjugate(x[i*incX], conjXt) * conjugate(y[i*incY], conjY);
    }
}

} } // namespace ref, ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_KERNEL_REF_DOTAXPY_TCC 1
