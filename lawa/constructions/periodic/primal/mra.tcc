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

namespace lawa {

template <typename T>
MRA<T,Primal,Periodic,CDF>::MRA(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j)
    : d(_d), d_(_d_), j0(j), mu(d&1), _j(j), phi(*this), M0(phi)
{

    // Refinement coefficients only in double prec. available due to missing support of higher
    // precision in blas routines.
    if (d==2) {

		// inner part
		_periodicRefCoeffs = new flens::DenseVector<flens::Array<long double> >[1];
		//_periodicRefCoeffs[0].engine().resize(3,0);
		//_periodicRefCoeffs[0] =  1.L/(2.L*std::sqrt(2.L)), 1.L/std::sqrt(2.L), 1.L/(2.L*std::sqrt(2.L));
		_periodicRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)1,0);
		_periodicRefCoeffs[0] =  1.L;
		_periodicOffsets = new FLENS_DEFAULT_INDEXTYPE[1];
        _periodicOffsets[0] =  1;

        _rightRefCoeffs = new flens::DenseVector<flens::Array<long double> >[1];
		_rightRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)2,0);
		_rightRefCoeffs[0] =  1.L, 1.L;
		_split = new FLENS_DEFAULT_INDEXTYPE[1];
		_split[0] = 1;
		_rightOffsets = new FLENS_DEFAULT_INDEXTYPE[1];
        _rightOffsets[0] =  1;

		_periodicL2Norms = new long double[1];
		_periodicL2Norms[0] =  std::sqrt(2.L/3.L);;
		_periodicH1SemiNorms = new long double[1];
		_periodicH1SemiNorms[0] = std::sqrt(2.L);

	}


}

template <typename T>
MRA<T,Primal,Periodic,CDF>::~MRA()
{
	if(d==2){
		delete[] _periodicRefCoeffs;
		delete[] _periodicOffsets;
		delete[] _rightRefCoeffs;
		delete[] _split;
		delete[] _rightOffsets;
		delete[] _periodicL2Norms;
		delete[] _periodicH1SemiNorms;
	}

}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,Periodic,CDF>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=0);
    _j = j;
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Periodic,CDF>::cardI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);   
    return pow2i<T>(j);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Periodic,CDF>::cardIL(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return std::max( std::ceil( (d - mu)/2.0 - 1), 0.0);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Periodic,CDF>::cardII(FLENS_DEFAULT_INDEXTYPE j) const
{
	assert(j >= j0);
    return pow2i<T>(j) - cardIL() - cardIR();
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,Periodic,CDF>::cardIR(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return std::ceil((d+mu)/2.0 - 1) + 1;
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Periodic,CDF>::rangeI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1,pow2i<T>(j));
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Periodic,CDF>::rangeIL(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    // If left index set is empty, return range = {1,0}
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1, cardIL());
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Periodic,CDF>::rangeII(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1 + cardIL() , pow2i<T>(j) - cardIR());
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Primal,Periodic,CDF>::rangeIR(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(pow2i<T>(j) - std::ceil((d+mu)/2.0 - 1), pow2i<T>(j));
}


} // namespace lawa

