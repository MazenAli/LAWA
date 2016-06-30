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

#include <lawa/methods/adaptive/datastructures/index.h>
#include <cassert>

namespace lawa {

Index1D::Index1D(void)
: j(0), k(0), xtype(XBSpline)
{}

//Index1D::Index1D(FLENS_DEFAULT_INDEXTYPE _j, FLENS_DEFAULT_INDEXTYPE _k, XType _xtype)
Index1D::Index1D(FLENS_DEFAULT_INDEXTYPE _j, FLENS_DEFAULT_INDEXTYPE _k, XType _xtype)
: j(_j), k(_k), xtype(_xtype)
{}

Index1D::Index1D(const Index1D &index)
: j(index.j), k(index.k), xtype(index.xtype)
{}

short
Index1D::levelSum(void) const
{
    return this->j;
}


Index1DC::Index1DC(void):
    Index1D(),
    d(0){}

Index1DC::Index1DC(FLENS_DEFAULT_INDEXTYPE _j, FLENS_DEFAULT_INDEXTYPE _k, XType _xtype, unsigned _d):
    Index1D(_j, _k, _xtype),
    d(_d){}

Index1DC::Index1DC(const Index1DC &index):
    Index1D(index),
    d(index.d){}


std::ostream& operator<<(std::ostream &s, const Index1D &_i)
{
    if (_i.xtype==XBSpline) {
        s << "scaling," << _i.j << "," << _i.k;
    } else {
        s << "wavelet," << _i.j << "," << _i.k;
    }
    return s;
}


std::ostream& operator<<(std::ostream &s, const Index1DC &_i)
{
    if (_i.xtype==XBSpline) {
        s << "scaling," << _i.j << "," << _i.k << "," << _i.d;
    } else {
        s << "wavelet," << _i.j << "," << _i.k << "," << _i.d;
    }
    return s;
}


Index2D::Index2D(void)
: index1(), index2()
{
}

Index2D::Index2D(const Index1D &_index1, const Index1D &_index2)
: index1(_index1), index2(_index2)
{
}

short
Index2D::levelSum(void) const
{
    return this->index1.j+this->index2.j;
}

std::ostream& operator<<(std::ostream &s, const Index2D &_i)
{
    s << _i.index1 << "," << _i.index2;
    return s;
}


Index3D::Index3D(void)
: index1(), index2(), index3()
{
}

Index3D::Index3D(const Index1D &_index1, const Index1D &_index2, const Index1D &_index3)
: index1(_index1), index2(_index2), index3(_index3)
{
}

Index3D::~Index3D(void)
{
}

short
Index3D::levelSum(void) const
{
    return this->index1.j+this->index2.j+this->index3.j;
}

std::ostream& operator<<(std::ostream &s, const Index3D &_i)
{
    s <<  _i.index1 << "," << _i.index2 << "," << _i.index3;
    return s;
}


IndexD::IndexD(void):
    dim_(0),
    index_(){}


IndexD::IndexD(const std::vector<Index1D>& _index):
    dim_(_index.size()),
    index_(_index){}


IndexD::IndexD(const size_type _dim):
    dim_(_dim),
    index_(_dim){}


FLENS_DEFAULT_INDEXTYPE
IndexD::levelSum() const
{
    FLENS_DEFAULT_INDEXTYPE ret = 0;
    for (size_type i=0; i<dim(); ++i) {
        ret += getIndex()[i].j;
    }

    return ret;
}


IndexD::size_type
IndexD::dim() const
{
    return dim_;
}


void
IndexD::setIndex(const size_type i, const Index1D& _index)
{
    assert(i>=1 && i<=dim());
    getIndex(i) = _index;
}


void
IndexD::setIndex(const std::vector<Index1D>& _index)
{
    dim_   = _index.size();
    index_ = _index;
}


const Index1D&
IndexD::getIndex(const size_type i) const
{
    assert(i>=1 && i<=dim());
    return index_[i-1];
}


Index1D&
IndexD::getIndex(const size_type i)
{
    assert(i>=1 && i<=dim());
    return index_[i-1];
}


const std::vector<Index1D>&
IndexD::getIndex() const
{
    return index_;
}


const Index1D&
IndexD::operator()(const size_type i) const
{
    assert(i>=1 && i<=dim());
    return getIndex(i);
}


Index1D&
IndexD::operator()(const size_type i)
{
    assert(i>=1 && i<=dim());
    return getIndex(i);
}


std::ostream& operator<<(std::ostream& s, const IndexD& _Index)
{
    typedef IndexD::size_type   size_type;
    s << "[";
    for (size_type i=1; i<=_Index.dim(); ++i) {
        s << _Index(i) << ";";
    }
    s << "]";

    return s;
}


} //namespace lawa
