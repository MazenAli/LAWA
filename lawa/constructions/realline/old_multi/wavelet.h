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

#ifndef LAWA_CONSTRUCTIONS_REALLINE_MULTI_WAVELET_H
#define LAWA_CONSTRUCTIONS_REALLINE_MULTI_WAVELET_H 1

#include <lawa/constructions/basisfunction.h>
#include <lawa/constructions/support.h>
#include <lawa/constructions/wavelet.h>

namespace lawa {

template <typename _T>
class Wavelet<_T,Orthogonal,R,Multi>
    : public BasisFunction<_T,Orthogonal,R,Multi>
{
    public:
        typedef _T T;
        static const FunctionSide Side = Orthogonal;
        static const DomainType Domain = R;
        static const Construction Cons = Multi;
    
        Wavelet(FLENS_DEFAULT_INDEXTYPE _d);
        
        Wavelet(const Basis<T,Orthogonal,R,Multi> &basis);
        
        virtual
        ~Wavelet();
        
        T
        operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const;
        
        Support<T>
        support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;
        
        Support<T>
        max_support() const;

        flens::DenseVector<flens::Array<T> >
        singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const;
        
        T
        tic(FLENS_DEFAULT_INDEXTYPE j) const;

//TODO        FLENS_DEFAULT_INDEXTYPE polynomialOrder;
        const FLENS_DEFAULT_INDEXTYPE d;
        const FLENS_DEFAULT_INDEXTYPE vanishingMoments;
        unsigned FLENS_DEFAULT_INDEXTYPE _numSplines;

    private:

        typedef T (*Evaluator)(T x, unsigned short deriv);
        
        FLENS_DEFAULT_INDEXTYPE
        _shift(FLENS_DEFAULT_INDEXTYPE k) const;
        
        FLENS_DEFAULT_INDEXTYPE
        _type(FLENS_DEFAULT_INDEXTYPE k) const;
    
        Evaluator *_evaluator;
        Support<T> *_support;
        flens::DenseVector<flens::Array<T> > *_singularSupport;

        Support<T> _max_support;

    
};
    
} // namespace lawa

#include <lawa/constructions/realline/multi/wavelet.tcc>

#endif // LAWA_CONSTRUCTIONS_REALLINE_MULTI_WAVELET_H
