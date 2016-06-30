/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

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
#include <lawa/constructions/interval/sparsemulti/_sparsemulti_scaling_evaluator.h>

namespace lawa {

template <typename T>
MRA<T,Primal,Interval,SparseMulti>::MRA(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE j)
    : d(_d), j0((j==-1) ? 0 : j), _bc(2,0), _j(j0), _numSplines(4), phi(*this)
{
    if (d!=4) {
        std::cerr << "MRA<T,Primal,Interval,SparseMulti> not yet implemented for d=" << d << std::endl;
        exit(1);
    }
    if (d==4) {
        if (j0<1) {
            std::cerr << "MRA<T,Primal,Interval,SparseMulti> not yet implemented for j0=" << j0 << std::endl;
            exit(1);
        }
    }
    setLevel(_j);
}

template <typename T>
MRA<T,Primal,Interval,SparseMulti>::~MRA()
{
    delete[] _innerEvaluator;
    delete[] _innerSupport;
    delete[] _innerSingularSupport;
}

template <typename T>
Support<T>
MRA<T,Primal,Interval,SparseMulti>::max_support() const
{
    if (d==4) return Support<T>(0.,2.);
    else{ // Control may reach end of non-void function
        assert(d==4);
        return Support<T>(-0.,0.);
    }
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::cardI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return 2*pow2i<FLENS_DEFAULT_INDEXTYPE>(j);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::cardIL(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return _numLeftParts;
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::cardII(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return 2*pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-_numLeftParts-_numRightParts;
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::cardIR(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return _numRightParts;
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Interval,SparseMulti>::rangeI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(0,cardI(j)-1);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Interval,SparseMulti>::rangeIL(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(0,0);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Interval,SparseMulti>::rangeII(FLENS_DEFAULT_INDEXTYPE j) const
{
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1,cardI(j)-2);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Interval,SparseMulti>::rangeIR(FLENS_DEFAULT_INDEXTYPE j) const
{
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(cardI(j)-1,cardI(j)-1);
}


//For adaptive schemes, we may require "FLENS_DEFAULT_INDEXTYPE" as index type for local scaling function repr.
template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::long_cardI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return 2*pow2i<FLENS_DEFAULT_INDEXTYPE>(j);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::long_cardIL(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return _numLeftParts;
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::long_cardII(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return 2*pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-_numLeftParts-_numRightParts;
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::long_cardIR(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return _numRightParts;
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Interval,SparseMulti>::long_rangeI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(0,long_cardI(j)-1);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Interval,SparseMulti>::long_rangeIL(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(0,0);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Interval,SparseMulti>::long_rangeII(FLENS_DEFAULT_INDEXTYPE j) const
{
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1,long_cardI(j)-2);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Interval,SparseMulti>::long_rangeIR(FLENS_DEFAULT_INDEXTYPE j) const
{
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(long_cardI(j)-1,long_cardI(j)-1);
}




template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Interval,SparseMulti>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,Interval,SparseMulti>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Primal,Interval,SparseMulti>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    switch (d) {
        case 4:
            _numSplines = 2;     //required for lambdaTilde
            // left B-splines
            _numLeftParts = 1;
            _leftEvaluator = new Evaluator[1];
            _leftEvaluator[0] = _cubic_sparsemulti_scaling_left_evaluator0;

            _leftSupport = new Support<T>[1];
            _leftSupport[0] = Support<T>(0,1);

            _leftSingularSupport = new flens::DenseVector<flens::Array<T> >[1];
            _leftSingularSupport[0] = linspace(0.,1.,2);

            _leftScalingFactors.engine().resize(1,0);
            _leftScalingFactors = 1.*std::sqrt(105.);

            // inner B-splines
            _numInnerParts = 2;
            _innerEvaluator = new Evaluator[2];
            _innerEvaluator[0] = _cubic_sparsemulti_scaling_inner_evaluator0;
            _innerEvaluator[1] = _cubic_sparsemulti_scaling_inner_evaluator1;

            _innerSupport = new Support<T>[2];
            _innerSupport[0] = Support<T>(-1,1);
            _innerSupport[1] = Support<T>(-1,1);

            _innerSingularSupport = new flens::DenseVector<flens::Array<T> >[2];
            _innerSingularSupport[0] = linspace(-1.,1.,3);
            _innerSingularSupport[1] = linspace(-1.,1.,3);

            _innerScalingFactors.engine().resize(2,0);
            _innerScalingFactors = std::sqrt(35./26.),std::sqrt(105./2.);

            // right B-splines
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _cubic_sparsemulti_scaling_right_evaluator0;

            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(-1,0);

            _rightSingularSupport = new flens::DenseVector<flens::Array<T> >[1];
            _rightSingularSupport[0] = linspace(-1.,0.,2);

            _rightScalingFactors.engine().resize(1,0);
            _rightScalingFactors = 1.*std::sqrt(105.);

            break;

        default: std::cerr << "BSpline<T,Primal,Interval,SparseMulti> not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }



    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;

}

} // namespace lawa
