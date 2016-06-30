/*
 LAWA - Library for Adaptive Wavelet Applications.
 Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <cassert>
#include <lawa/math/binomial.h>

namespace lawa {

FLENS_DEFAULT_INDEXTYPE
binomial(FLENS_DEFAULT_INDEXTYPE n, FLENS_DEFAULT_INDEXTYPE k)
{
    if (k>n) {
        return 0;
    }
    if (k==n) {
        return 1;
    }
    if (k==0) {
        return 1;
    }
    FLENS_DEFAULT_INDEXTYPE res = 1;
    for (FLENS_DEFAULT_INDEXTYPE i=k+1; i<=n; ++i) {
        res *= i;
    }
    FLENS_DEFAULT_INDEXTYPE den = 1;
    for (FLENS_DEFAULT_INDEXTYPE i=2; i<=n-k; ++i) {
        den *= i;
    }
    assert((1.*res)/den == res/den);
    return res / den;
}

} // namespace lawa

