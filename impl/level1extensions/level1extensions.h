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

#ifndef ULMBLAS_IMPL_LEVEL1EXTENSIONS_LEVEL1EXTENSIONS_H
#define ULMBLAS_IMPL_LEVEL1EXTENSIONS_LEVEL1EXTENSIONS_H 1

#include <ulmblas/impl/level1extensions/asum1.h>
#include <ulmblas/impl/level1extensions/axpby.h>
#include <ulmblas/impl/level1extensions/axpy2v.h>
#include <ulmblas/impl/level1extensions/axpyf.h>
#include <ulmblas/impl/level1extensions/conj.h>
#include <ulmblas/impl/level1extensions/dot2v.h>
#include <ulmblas/impl/level1extensions/dotaxpy.h>
#include <ulmblas/impl/level1extensions/dotxaxpyf.h>
#include <ulmblas/impl/level1extensions/dotxf.h>
#include <ulmblas/impl/level1extensions/geaxpy.h>
#include <ulmblas/impl/level1extensions/gecopy.h>
#include <ulmblas/impl/level1extensions/geconj.h>
#include <ulmblas/impl/level1extensions/gecotr.h>
#include <ulmblas/impl/level1extensions/geraxpy.h>
#include <ulmblas/impl/level1extensions/gescal.h>
#include <ulmblas/impl/level1extensions/geswap.h>
#include <ulmblas/impl/level1extensions/iamax1.h>
#include <ulmblas/impl/level1extensions/kernel/axpy2v.h>
#include <ulmblas/impl/level1extensions/kernel/axpyf.h>
#include <ulmblas/impl/level1extensions/kernel/dot2v.h>
#include <ulmblas/impl/level1extensions/kernel/dotaxpy.h>
#include <ulmblas/impl/level1extensions/kernel/dotxaxpyf.h>
#include <ulmblas/impl/level1extensions/kernel/dotxf.h>
#include <ulmblas/impl/level1extensions/raxpy.h>
#include <ulmblas/impl/level1extensions/trlaxpy.h>
#include <ulmblas/impl/level1extensions/trlcopy.h>
#include <ulmblas/impl/level1extensions/trlscal.h>
#include <ulmblas/impl/level1extensions/truaxpy.h>
#include <ulmblas/impl/level1extensions/truscal.h>

#endif // ULMBLAS_IMPL_LEVEL1EXTENSIONS_LEVEL1EXTENSIONS_H
