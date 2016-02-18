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

#ifndef ULMBLAS_IMPL_AUXILIARY_MEMORYPOOL_TCC
#define ULMBLAS_IMPL_AUXILIARY_MEMORYPOOL_TCC 1

#include <cassert>
#include <list>
#include <unordered_map>
#include <type_traits>

#include <ulmblas/impl/auxiliary/isfundamental.h>
#include <ulmblas/impl/auxiliary/memorypool.h>
#include <ulmblas/impl/auxiliary/malloc_aligned.h>
#include <ulmblas/impl/config/blocksize.h>


namespace ulmBLAS {

template <typename T>
T *
MemoryPool<T>::allocate(size_t n)
{
    mutex_.lock();

    constexpr std::size_t align = BlockSize<T>::align;

    BlockList &free = free_[n];
    T         *block;

    if (free.empty()) {
        if (align!=0) {
            block = malloc_aligned<T,align>(n);
        } else {
            block = new T[n];
        }
        allocated_.push_back(block);
    } else {
        block = free.back();
        free.pop_back();
    }
    used_[block] = n;

    mutex_.unlock();
    return block;
}

template <typename T>
void
MemoryPool<T>::release(T *block)
{
    mutex_.lock();

    if (block) {
        assert(used_.count(block)==1);
        size_t n = used_[block];
        free_[n].push_back(block);
    }

    mutex_.unlock();
}

template <typename T>
MemoryPool<T>::~MemoryPool()
{
    constexpr std::size_t align = BlockSize<T>::align;

    while (!allocated_.empty()) {
        if (align!=0) {
            free_aligned(allocated_.back());
        } else {
            delete [] allocated_.back();
        }
        allocated_.pop_back();
    }
}

} // namespace ulmBLAS

#endif // ULMBLAS_IMPL_AUXILIARY_MEMORYPOOL_TCC
