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

#ifndef ULMBLAS_IMPL_CONFIG_BLOCKSIZE_H
#define ULMBLAS_IMPL_CONFIG_BLOCKSIZE_H 1

#include <complex>
#include <type_traits>

namespace ulmBLAS {

#if defined(USE_TESTPARAM)

#   define BS_DEFAULT_MC_S          USE_TESTPARAM_MC
#   define BS_DEFAULT_KC_S          USE_TESTPARAM_KC
#   define BS_DEFAULT_NC_S          USE_TESTPARAM_NC
#   define BS_DEFAULT_MR_S          USE_TESTPARAM_MR
#   define BS_DEFAULT_NR_S          USE_TESTPARAM_NR

#   define BS_DEFAULT_MC_D          USE_TESTPARAM_MC
#   define BS_DEFAULT_KC_D          USE_TESTPARAM_KC
#   define BS_DEFAULT_NC_D          USE_TESTPARAM_NC
#   define BS_DEFAULT_MR_D          USE_TESTPARAM_MR
#   define BS_DEFAULT_NR_D          USE_TESTPARAM_NR

#   define BS_DEFAULT_MC_C          USE_TESTPARAM_MC
#   define BS_DEFAULT_KC_C          USE_TESTPARAM_KC
#   define BS_DEFAULT_NC_C          USE_TESTPARAM_NC
#   define BS_DEFAULT_MR_C          USE_TESTPARAM_MR
#   define BS_DEFAULT_NR_C          USE_TESTPARAM_NR

#   define BS_DEFAULT_MC_Z          USE_TESTPARAM_MC
#   define BS_DEFAULT_KC_Z          USE_TESTPARAM_KC
#   define BS_DEFAULT_NC_Z          USE_TESTPARAM_NC
#   define BS_DEFAULT_MR_Z          USE_TESTPARAM_MR
#   define BS_DEFAULT_NR_Z          USE_TESTPARAM_NR

#elif defined(USE_FMA)

#   define SIMD_WIDTH               256

#   define BS_DEFAULT_MC_S          144
#   define BS_DEFAULT_KC_S          256
#   define BS_DEFAULT_NC_S          4080
#   define BS_DEFAULT_MR_S          16
#   define BS_DEFAULT_NR_S          6

#   ifdef USE_4_12
#      define BS_DEFAULT_MC_D       72
#      define BS_DEFAULT_KC_D       256
#      define BS_DEFAULT_NC_D       4080
#      define BS_DEFAULT_MR_D       8
#      define BS_DEFAULT_NR_D       6
#   else
#      define BS_DEFAULT_MC_D       96
#      define BS_DEFAULT_KC_D       192
#      define BS_DEFAULT_NC_D       4080
#      define BS_DEFAULT_MR_D       4
#      define BS_DEFAULT_NR_D       12
#   endif

#   define BS_DEFAULT_MC_C          96
#   define BS_DEFAULT_KC_C          256
#   define BS_DEFAULT_NC_C          4096
#   define BS_DEFAULT_MR_C          8
#   define BS_DEFAULT_NR_C          4

#   define BS_DEFAULT_MC_Z          64
#   define BS_DEFAULT_KC_Z          192
#   define BS_DEFAULT_NC_Z          4096
#   define BS_DEFAULT_MR_Z          4
#   define BS_DEFAULT_NR_Z          4


#elif defined(USE_AVX)

#   define SIMD_WIDTH               256

#   define BS_DEFAULT_MC_S          128
#   define BS_DEFAULT_KC_S          384
#   define BS_DEFAULT_NC_S          4096
#   define BS_DEFAULT_MR_S          8
#   define BS_DEFAULT_NR_S          8

#   define BS_DEFAULT_MC_D          96
#   define BS_DEFAULT_KC_D          256
#   define BS_DEFAULT_NC_D          4096
#   define BS_DEFAULT_MR_D          4
#   define BS_DEFAULT_NR_D          8

#   define BS_DEFAULT_MC_C          96
#   define BS_DEFAULT_KC_C          256
#   define BS_DEFAULT_NC_C          4096
#   define BS_DEFAULT_MR_C          8
#   define BS_DEFAULT_NR_C          4

#   define BS_DEFAULT_MC_Z          64
#   define BS_DEFAULT_KC_Z          192
#   define BS_DEFAULT_NC_Z          4096
#   define BS_DEFAULT_MR_Z          4
#   define BS_DEFAULT_NR_Z          4

#elif defined(USE_SSE)

#   define SIMD_WIDTH               128

#   define BS_DEFAULT_MC_S          768
#   define BS_DEFAULT_KC_S          384
#   define BS_DEFAULT_NC_S          4096
#   define BS_DEFAULT_MR_S          8
#   define BS_DEFAULT_NR_S          4

#   define BS_DEFAULT_MC_D          384
#   define BS_DEFAULT_KC_D          384
#   define BS_DEFAULT_NC_D          4096
#   define BS_DEFAULT_MR_D          4
#   define BS_DEFAULT_NR_D          4

#   define BS_DEFAULT_MC_C          384
#   define BS_DEFAULT_KC_C          384
#   define BS_DEFAULT_NC_C          4096
#   define BS_DEFAULT_MR_C          4
#   define BS_DEFAULT_NR_C          2

#   define BS_DEFAULT_MC_Z          192
#   define BS_DEFAULT_KC_Z          384
#   define BS_DEFAULT_NC_Z          4096
#   define BS_DEFAULT_MR_Z          2
#   define BS_DEFAULT_NR_Z          2

#else

#   define SIMD_WIDTH               0

#   define BS_DEFAULT_MC_S          768
#   define BS_DEFAULT_KC_S          384
#   define BS_DEFAULT_NC_S          4096
#   define BS_DEFAULT_MR_S          8
#   define BS_DEFAULT_NR_S          4

#   define BS_DEFAULT_MC_D          384
#   define BS_DEFAULT_KC_D          384
#   define BS_DEFAULT_NC_D          4096
#   define BS_DEFAULT_MR_D          4
#   define BS_DEFAULT_NR_D          4

#   define BS_DEFAULT_MC_C          384
#   define BS_DEFAULT_KC_C          384
#   define BS_DEFAULT_NC_C          4096
#   define BS_DEFAULT_MR_C          4
#   define BS_DEFAULT_NR_C          2

#   define BS_DEFAULT_MC_Z          192
#   define BS_DEFAULT_KC_Z          384
#   define BS_DEFAULT_NC_Z          4096
#   define BS_DEFAULT_MR_Z          2
#   define BS_DEFAULT_NR_Z          2

#endif

template <typename T>
struct BlockSize
{
    typedef float                   S;
    typedef double                  D;
    typedef std::complex<float>     C;
    typedef std::complex<double>    Z;

    static const int MC = std::is_same<S, T>::value ? BS_DEFAULT_MC_S
                        : std::is_same<D, T>::value ? BS_DEFAULT_MC_D
                        : std::is_same<C, T>::value ? BS_DEFAULT_MC_C
                        : std::is_same<Z, T>::value ? BS_DEFAULT_MC_Z
                        : 192;
    static const int KC = std::is_same<S, T>::value ? BS_DEFAULT_KC_S
                        : std::is_same<D, T>::value ? BS_DEFAULT_KC_D
                        : std::is_same<C, T>::value ? BS_DEFAULT_KC_C
                        : std::is_same<Z, T>::value ? BS_DEFAULT_KC_Z
                        : 192;
    static const int NC = std::is_same<S, T>::value ? BS_DEFAULT_NC_S
                        : std::is_same<D, T>::value ? BS_DEFAULT_NC_D
                        : std::is_same<C, T>::value ? BS_DEFAULT_NC_C
                        : std::is_same<Z, T>::value ? BS_DEFAULT_NC_Z
                        : 4096;
    static const int MR = std::is_same<S, T>::value ? BS_DEFAULT_MR_S
                        : std::is_same<D, T>::value ? BS_DEFAULT_MR_D
                        : std::is_same<C, T>::value ? BS_DEFAULT_MR_C
                        : std::is_same<D, T>::value ? BS_DEFAULT_MR_Z
                        : 2;

    static const int NR = std::is_same<S, T>::value ? BS_DEFAULT_NR_S
                        : std::is_same<D, T>::value ? BS_DEFAULT_NR_D
                        : std::is_same<C, T>::value ? BS_DEFAULT_NR_C
                        : std::is_same<Z, T>::value ? BS_DEFAULT_NR_Z
                        : 2;

    static constexpr int rwidth = SIMD_WIDTH;   // SIMD-Register width in bits
    static constexpr int align  = rwidth / 8;   // SIMD-Register width in bytes
    static constexpr int vlen   = rwidth / (8*sizeof(T));

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
};



} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_CONFIG_BLOCKSIZE_H 1
