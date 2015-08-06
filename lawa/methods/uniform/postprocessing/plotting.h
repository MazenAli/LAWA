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
 
#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_PLOTTING_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_PLOTTING_H 1

namespace lawa {

template<typename T, typename Basis1D>
void 
print_U(const flens::DenseVector<flens::Array<T> >& u, const Basis1D& basis, const int J, 
        const char* filename, const double deltaX=1./128.);
        
template<typename T, typename Basis2D>
void
print_U(const flens::DenseVector<flens::Array<T> >& u, const Basis2D& basis, const int J_x, const int J_y, 
                   const char* filename, const double deltaX=0.01, const double deltaY=0.01);
        
template<typename T, typename Basis1D>
void 
print_U(const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U, const Basis1D& basis, 
        const int J, const char* filename, const T timestep, const int K, const double deltaX=1./128.);
    
} // namespace lawa

#include <lawa/methods/uniform/postprocessing/plotting.tcc>

#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_PLOTTING_H
