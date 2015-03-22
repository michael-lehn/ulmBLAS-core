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

#ifndef ULMBLAS_IMPL_LEVEL1_SWAP_TCC
#define ULMBLAS_IMPL_LEVEL1_SWAP_TCC 1

#include <utility>
#include <ulmblas/impl/level1/swap.h>
#include <ulmblas/impl/auxiliary/conjugate.h>

namespace ulmBLAS {

template <typename IndexType, typename VX, typename VY>
void
swap(IndexType      n,
     VX             *x,
     IndexType      incX,
     VY             *y,
     IndexType      incY)
{
    const IndexType    UnitStride(1);

    if (incX==UnitStride && incY==UnitStride) {
        for (IndexType i=0; i<n; ++i) {
            std::swap(x[i], y[i]);
        }
     } else {
        for (IndexType i=0; i<n; ++i) {
            std::swap(x[i*incX], y[i*incY]);
        }
    }
}

template <typename IndexType, typename VX, typename VY>
void
swapc(IndexType      n,
      VX             *x,
      IndexType      incX,
      VY             *y,
      IndexType      incY)
{
    const IndexType    UnitStride(1);

    if (incX==UnitStride && incY==UnitStride) {
        for (IndexType i=0; i<n; ++i) {
            std::swap(x[i], y[i]);
            x[i] = conjugate(x[i]);
            y[i] = conjugate(y[i]);
        }
     } else {
        for (IndexType i=0; i<n; ++i) {
            std::swap(x[i*incX], y[i*incY]);
            x[i*incX] = conjugate(x[i*incX]);
            y[i*incY] = conjugate(y[i*incY]);
        }
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_LEVEL1_SWAP_TCC 1
