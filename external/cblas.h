/*
 *   Copyright (c) 2015, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ULMBLAS_EXTERNAL_CBLAS_H
#define ULMBLAS_EXTERNAL_CBLAS_H 1

#include <cassert>

namespace cblas {

#ifdef WITH_MKL
#   define WITH_CBLAS
#   include <mkl_cblas.h>
#   define BLAS_IMPL   "Intel MKL"
    typedef MKL_INT  CBLAS_INT;
#endif

#ifdef WITH_CBLAS

typedef std::complex<float>     ComplexFloat;
typedef std::complex<double>    ComplexDouble;

#ifdef CBLAS_DEBUG
#   ifndef CBLAS_DEBUG_OUT
#   define CBLAS_DEBUG_OUT(msg)  std::cerr << msg << std::endl
#   endif
#else
#   ifndef CBLAS_DEBUG_OUT
#   define CBLAS_DEBUG_OUT(msg)
#   endif
#endif

template <typename VOID=void>
struct Blas
{
    static const CBLAS_ORDER
    order(bool colMajor)
    {
        if (colMajor) {
            return CblasColMajor;
        }
        return CblasRowMajor;
    }

    static const CBLAS_TRANSPOSE
    trans(bool trans, bool conj)
    {
        if (!trans && !conj) {
            return CblasNoTrans;
        } else if (trans && !conj) {
            return CblasTrans;
        } else if (trans && conj) {
            return CblasConjTrans;
        } else {
            assert(0);
            return CblasNoTrans;
        }
    }

    static const CBLAS_SIDE
    side(bool left)
    {
        return left ? CblasLeft : CblasRight;
    }

    static const CBLAS_UPLO
    uplo(bool lower)
    {
        return lower ? CblasLower : CblasUpper;
    }

    static const CBLAS_DIAG
    diag(bool unitDiag)
    {
        return unitDiag ? CblasUnit : CblasNonUnit;
    }
};

template <typename IndexType, typename ElementType=void, typename Enable=void>
struct If
{
};

template <typename IndexType>
struct If<IndexType, void,
          typename std::enable_if<
                        std::is_convertible<IndexType, CBLAS_INT>::value
                        >::type>
{
    typedef void isBlasCompatibleInteger;
};

template <typename IndexType, typename ElementType>
struct If<IndexType, ElementType,
          typename std::enable_if<
                        std::is_convertible<IndexType, CBLAS_INT>::value
                    && (std::is_same<ElementType, float>::value ||
                        std::is_same<ElementType, double>::value ||
                        std::is_same<ElementType, ComplexFloat>::value ||
                        std::is_same<ElementType, ComplexDouble>::value)
                        >::type>
{
    typedef void isBlasCompatible;
};

#endif // WITH_CBLAS

} // namespace cblas

#endif // ULMBLAS_EXTERNAL_CBLAS_H 1
