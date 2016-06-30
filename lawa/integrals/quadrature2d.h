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

#ifndef LAWA_INTEGRALS_QUADRATURE2D_H
#define LAWA_INTEGRALS_QUADRATURE2D_H 1

#include <extensions/sparsegrid/sparse_grid_mixed.h>
#include <lawa/math/lawa_math.h>
#include <lawa/settings/enum.h>
#include <lawa/integrals/integral2d.h>

namespace lawa {
  
template<QuadratureType, typename, typename> class Integral2D;

template <QuadratureType Quad, typename _Integral2D>
    class Quadrature2D;
//{
//};

template <typename _Integral2D>
class Quadrature2D<SparseGridGP,_Integral2D>
{
    public:
        //typedef typename _Integral2D::T T;
        typedef double T;

        Quadrature2D(const _Integral2D &integral);

        double
        operator()(double ax, double bx, double ay, double by) const;

        void setOrder(FLENS_DEFAULT_INDEXTYPE order);
        void setLevel(FLENS_DEFAULT_INDEXTYPE level);

        FLENS_DEFAULT_INDEXTYPE numGridPoints;

    private:
        const _Integral2D &_integral;
        FLENS_DEFAULT_INDEXTYPE _level;

        void
        _initSparseGrid();

        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;
};

template <typename _Integral2D>
class Quadrature2D<FullGridGL,_Integral2D>
{
    public:
        typedef typename _Integral2D::T T;

        Quadrature2D(const _Integral2D &integral);

        const T
        operator()(T ax, T bx, T ay, T by) const;

        void setOrder(FLENS_DEFAULT_INDEXTYPE order);

        void _initFullGrid();
        const _Integral2D &_integral;

    private:
        FLENS_DEFAULT_INDEXTYPE _order;

        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _knots;
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > _weights;
};

template <typename _Integral2D>
class Quadrature2D<FullGridGL_localOrder, _Integral2D>
{
    typedef Integral2D<FullGridGL, typename _Integral2D::Basis_X, typename _Integral2D::Basis_Y> Integral_nonLocal;

    public:
        typedef typename _Integral2D::T T;

        Quadrature2D(const _Integral2D &integral);

        const T
        operator()(T ax, T bx, T ay, T by) const;

        void setOrder(FLENS_DEFAULT_INDEXTYPE order);

        void set_refindfct(T (*_fct)(T));

        void set_refindtol(T _tol);

        void set_lowOrder(FLENS_DEFAULT_INDEXTYPE _lowOrder);

        void set_highOrder(FLENS_DEFAULT_INDEXTYPE _highOrder);

    private:

        const _Integral2D &_integral;
        Integral_nonLocal integral_nonLocal;

        Quadrature2D<FullGridGL, Integral_nonLocal> quadrature_lowOrder;
        Quadrature2D<FullGridGL, Integral_nonLocal> quadrature_highOrder;
        T (*refindfct)(T);
        T refindtol;
        FLENS_DEFAULT_INDEXTYPE lowOrder, highOrder;

        bool
        ref_indicator(T at, T bt, T ax, T bx) const;
};

}   //namespace lawa

#include <lawa/integrals/quadrature2d.tcc>

#endif   //LAWA_RIGHTHANDSIDES_QUADRATURE2D_H

