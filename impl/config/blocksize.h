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

#ifndef ULMBLAS_IMPL_CONFIG_BLOCKSIZE_H
#define ULMBLAS_IMPL_CONFIG_BLOCKSIZE_H 1

#include <complex>
#include <type_traits>

#include <ulmblas/impl/config/simd.h>

namespace ulmBLAS {

#if defined(USE_TESTPARAM)

template <typename T>
struct BlockSize
{
    static const int MC = USE_TESTPARAM_MC;
    static const int KC = USE_TESTPARAM_KC;
    static const int NC = USE_TESTPARAM_NC;


    static const int MR = USE_TESTPARAM_MR;
    static const int NR = USE_TESTPARAM_NR;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
};

#elif defined(USE_SSE)

template <typename T>
struct BlockSize
{
    static const int MC = std::is_same<float, T>::value                ? 768
                        : std::is_same<double, T>::value               ? 384
                        : std::is_same<std::complex<float>, T>::value  ? 384
                        : std::is_same<std::complex<double>, T>::value ? 192
                        : -1;

    static const int KC = std::is_same<float, T>::value                ? 384
                        : std::is_same<double, T>::value               ? 384
                        : std::is_same<std::complex<float>, T>::value  ? 384
                        : std::is_same<std::complex<double>, T>::value ? 192
                        : -1;

    static const int NC = std::is_same<float, T>::value                ? 4096
                        : std::is_same<double, T>::value               ? 4096
                        : std::is_same<std::complex<float>, T>::value  ? 4096
                        : std::is_same<std::complex<double>, T>::value ? 4096
                        : -1;


    static const int MR = std::is_same<float, T>::value                ? 8
                        : std::is_same<double, T>::value               ? 4
                        : std::is_same<std::complex<float>, T>::value  ? 4
                        : std::is_same<std::complex<double>, T>::value ? 2
                        : -1;

    static const int NR = std::is_same<float, T>::value                ? 4
                        : std::is_same<double, T>::value               ? 4
                        : std::is_same<std::complex<float>, T>::value  ? 2
                        : std::is_same<std::complex<double>, T>::value ? 2
                        : -1;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
};

#else

template <typename T>
struct BlockSize
{
    static const int MC = std::is_same<float, T>::value                ? 768
                        : std::is_same<double, T>::value               ? 384
                        : std::is_same<std::complex<float>, T>::value  ? 384
                        : std::is_same<std::complex<double>, T>::value ? 192
                        : -1;

    static const int KC = std::is_same<float, T>::value                ? 384
                        : std::is_same<double, T>::value               ? 384
                        : std::is_same<std::complex<float>, T>::value  ? 384
                        : std::is_same<std::complex<double>, T>::value ? 192
                        : -1;

    static const int NC = std::is_same<float, T>::value                ? 4096
                        : std::is_same<double, T>::value               ? 4096
                        : std::is_same<std::complex<float>, T>::value  ? 4096
                        : std::is_same<std::complex<double>, T>::value ? 4096
                        : -1;


    static const int MR = std::is_same<float, T>::value                ? 8
                        : std::is_same<double, T>::value               ? 4
                        : std::is_same<std::complex<float>, T>::value  ? 4
                        : std::is_same<std::complex<double>, T>::value ? 2
                        : -1;

    static const int NR = std::is_same<float, T>::value                ? 4
                        : std::is_same<double, T>::value               ? 4
                        : std::is_same<std::complex<float>, T>::value  ? 2
                        : std::is_same<std::complex<double>, T>::value ? 2
                        : -1;

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
};

#endif


} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_CONFIG_BLOCKSIZE_H 1
