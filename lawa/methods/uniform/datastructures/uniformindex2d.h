/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2014  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

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

#ifndef LAWA_METHODS_UNIFORM_DATASTRUCTURES_UNIFORMINDEX2D_H
#define LAWA_METHODS_UNIFORM_DATASTRUCTURES_UNIFORMINDEX2D_H 1

#include <lawa/settings/enum.h>

namespace lawa {
  
template<typename Basis>    
class UniformIndex2D
{
    private:
        const Basis& basis;
        const FLENS_DEFAULT_INDEXTYPE J_x, J_y;
        
        const FLENS_DEFAULT_INDEXTYPE offsetIx;
        const FLENS_DEFAULT_INDEXTYPE offsetIy;
        const FLENS_DEFAULT_INDEXTYPE offsetJx;
        const FLENS_DEFAULT_INDEXTYPE offsetJy;
    
    public:
        UniformIndex2D(const Basis& _basis, const FLENS_DEFAULT_INDEXTYPE _J_x, const FLENS_DEFAULT_INDEXTYPE _J_y);
        
        FLENS_DEFAULT_INDEXTYPE
        operator()(XType xtype_x, FLENS_DEFAULT_INDEXTYPE jx, FLENS_DEFAULT_INDEXTYPE kx,
                   XType xtype_y, FLENS_DEFAULT_INDEXTYPE jy, FLENS_DEFAULT_INDEXTYPE ky) const;
};
    
    
} // namespace lawa

#include <lawa/methods/uniform/datastructures/uniformindex2d.tcc>

#endif // LAWA_METHODS_UNIFORM_DATASTRUCTURES_UNIFORMINDEX2D_H

