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
 
#ifndef LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_BASIS_H
#define LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_BASIS_H 1

#include <cassert>
#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/bspline.h>
#include <lawa/constructions/mra.h>
#include <lawa/constructions/wavelet.h>
#include <lawa/constructions/interval/sparsemulti/_sparsemulti_wavelet_evaluator.h>


namespace lawa {

template <typename _T>
class Basis<_T,Primal,RPlus,SparseMulti>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Primal;
        static const DomainType Domain = RPlus;
        static const Construction Cons = SparseMulti;
        
        typedef BasisFunction<T,Primal,RPlus,SparseMulti> BasisFunctionType;
        typedef BSpline<T,Primal,RPlus,SparseMulti>       BSplineType;
        typedef Wavelet<T,Primal,RPlus,SparseMulti>       WaveletType;

        Basis(const int d, const int j=-1);
    
        virtual
        ~Basis();
    
        int
        level() const;
    
        void
        setLevel(const int j) const;
    
        template <BoundaryCondition BC>
            void
            enforceBoundaryCondition();
    
        const BasisFunctionType &
        generator(XType xtype) const;

        //--- cardinalities of left index set.
        long
        cardJL(const int j) const;

        //--- ranges of left index set.
        const flens::Range<long>
        rangeJL(const int j) const;
    
        MRA<T,Primal,RPlus,SparseMulti> mra;
    
        const int d;
        const int j0;          // minimal used(!) level.
    
    private:
        flens::DenseVector<flens::Array<int> > _bc;  // the boundary conditions
                                       // bc(0) = 1 -> Dirichlet BC left.
        
        mutable int _j;                // the current level.
    
        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        friend class Wavelet<T,Primal,RPlus,SparseMulti>;

        unsigned int _numLeftParts,
                     _numInnerParts;
        Evaluator *_leftEvaluator, 
                  *_innerEvaluator;
        Support<T> *_leftSupport, 
                   *_innerSupport;
        flens::DenseVector<flens::Array<T> > *_leftSingularSupport, 
                               *_innerSingularSupport;
        flens::DenseVector<flens::Array<T> > _leftScalingFactors, _innerScalingFactors;
        
    public:
        Wavelet<T,Primal,RPlus,SparseMulti> psi;
};

} // namespace lawa

#include <lawa/constructions/rplus/sparsemulti/basis.tcc>

#endif // LAWA_CONSTRUCTIONS_RPLUS_SPARSEMULTI_BASIS_H
