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

#ifndef ULMBLAS_IMPL_LEVEL2_LEVEL2_H
#define ULMBLAS_IMPL_LEVEL2_LEVEL2_H 1

#include <ulmblas/impl/level2/gbmtv.h>
#include <ulmblas/impl/level2/gbmv.h>
#include <ulmblas/impl/level2/gemv.h>
#include <ulmblas/impl/level2/ger.h>
#include <ulmblas/impl/level2/hblmv.h>
#include <ulmblas/impl/level2/hbumv.h>
#include <ulmblas/impl/level2/helmv.h>
#include <ulmblas/impl/level2/helr.h>
#include <ulmblas/impl/level2/helr2.h>
#include <ulmblas/impl/level2/hplmv.h>
#include <ulmblas/impl/level2/hplr.h>
#include <ulmblas/impl/level2/hplr2.h>
#include <ulmblas/impl/level2/hpumv.h>
#include <ulmblas/impl/level2/hpur.h>
#include <ulmblas/impl/level2/hpur2.h>
#include <ulmblas/impl/level2/level2.h>
#include <ulmblas/impl/level2/sblmv.h>
#include <ulmblas/impl/level2/sbumv.h>
#include <ulmblas/impl/level2/splmv.h>
#include <ulmblas/impl/level2/splr.h>
#include <ulmblas/impl/level2/splr2.h>
#include <ulmblas/impl/level2/spumv.h>
#include <ulmblas/impl/level2/spur.h>
#include <ulmblas/impl/level2/spur2.h>
#include <ulmblas/impl/level2/sylmv.h>
#include <ulmblas/impl/level2/sylr.h>
#include <ulmblas/impl/level2/sylr2.h>
#include <ulmblas/impl/level2/tblmtv.h>
#include <ulmblas/impl/level2/tblmv.h>
#include <ulmblas/impl/level2/tblstv.h>
#include <ulmblas/impl/level2/tblsv.h>
#include <ulmblas/impl/level2/tbumtv.h>
#include <ulmblas/impl/level2/tbumv.h>
#include <ulmblas/impl/level2/tbustv.h>
#include <ulmblas/impl/level2/tbusv.h>
#include <ulmblas/impl/level2/tplmtv.h>
#include <ulmblas/impl/level2/tplmv.h>
#include <ulmblas/impl/level2/tplstv.h>
#include <ulmblas/impl/level2/tplsv.h>
#include <ulmblas/impl/level2/tpumtv.h>
#include <ulmblas/impl/level2/tpumv.h>
#include <ulmblas/impl/level2/tpustv.h>
#include <ulmblas/impl/level2/tpusv.h>
#include <ulmblas/impl/level2/trlmv.h>
#include <ulmblas/impl/level2/trlsv.h>
#include <ulmblas/impl/level2/trumv.h>
#include <ulmblas/impl/level2/trusv.h>

#endif // ULMBLAS_IMPL_LEVEL2_LEVEL2_H
