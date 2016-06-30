/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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
 
namespace lawa{
    
template<typename FirstBasis, typename SecondBasis>
TensorBasis2D<Uniform, FirstBasis, SecondBasis>::TensorBasis2D(const FirstBasis &_basis1, const SecondBasis &_basis2)  
    : first(_basis1), second(_basis2)
{
}

template<typename FirstBasis, typename SecondBasis>
FLENS_DEFAULT_INDEXTYPE
TensorBasis2D<Uniform, FirstBasis, SecondBasis>::dim(FLENS_DEFAULT_INDEXTYPE J_x, FLENS_DEFAULT_INDEXTYPE J_y) const
{
    return first.mra.cardI(J_x) * second.mra.cardI(J_y);
}

template<typename FirstBasis, typename SecondBasis>
FLENS_DEFAULT_INDEXTYPE
TensorBasis2D<Uniform, FirstBasis, SecondBasis>::J1_max(const FLENS_DEFAULT_INDEXTYPE J_x, const FLENS_DEFAULT_INDEXTYPE /*J_y*/, const FLENS_DEFAULT_INDEXTYPE /*jy*/) const
{
    return J_x;
}

template<typename FirstBasis, typename SecondBasis>
FLENS_DEFAULT_INDEXTYPE
TensorBasis2D<Uniform, FirstBasis, SecondBasis>::J2_max(const FLENS_DEFAULT_INDEXTYPE /*J_x*/, const FLENS_DEFAULT_INDEXTYPE J_y, const FLENS_DEFAULT_INDEXTYPE /*jx*/) const
{
    return J_y;
}

} // namespace lawa
