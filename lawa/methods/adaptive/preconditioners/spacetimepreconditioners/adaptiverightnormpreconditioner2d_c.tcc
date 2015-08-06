#include <cmath>

namespace lawa {

template <typename T, typename Basis2D>
AdaptiveRightNormPreconditioner2D_c<T,Basis2D>::AdaptiveRightNormPreconditioner2D_c(const Basis2D &basis, T s)
    : _s(s), _integral_t(basis.first,basis.first), 
             _integral_x(basis.second,basis.second)
{
}

/* Do _all_ calculations here (and only redirect from the operator()(j1,k1,xtype1,..):
 * 		Usually, we call the preconditioner in an adaptive setting using the Index2D, so
 * 		that this saves one function call.
 */
template <typename T, typename Basis2D>
T
AdaptiveRightNormPreconditioner2D_c<T,Basis2D>::operator()(const Index2D &index)
{
	T value_t, value_x;

	// See if we have calculated the entry in time before
	auto it_t = values_L2_t.find(index.index1);
	// If yes, reuse the values
	if(it_t != values_L2_t.end()){
		value_t = (*it_t).second;
	}
	else{
		value_t = _integral_t(index.index1.j,index.index1.k,index.index1.xtype,0,
							  index.index1.j,index.index1.k,index.index1.xtype,0);
		values_L2_t.insert(std::make_pair(index.index1,value_t));
	}

	// See if we have calculated the entry in space before
	auto it_x = values_L2_x.find(index.index2);
	// If yes, reuse the values
	if(it_x != values_L2_x.end()){
		value_x = (*it_x).second;
	}
	else{
		value_x = _integral_x(index.index2.j,index.index2.k,index.index2.xtype,0,
							  index.index2.j,index.index2.k,index.index2.xtype,0);
		values_L2_x.insert(std::make_pair(index.index2,value_x));
	}

	if(_s==2){
		auto it_H1 = values_H1semi.find(index.index2);
		if(it_H1 != values_H1semi.end()){
			return 1./std::sqrt( value_t*(value_x+(*it_H1).second) + pow2i<T>(2*index.index1.j)*value_t/(value_x + (*it_H1).second));
		}
		else{
			T dd_value_x = _integral_x(index.index2.j,index.index2.k,index.index2.xtype,1,
					  	  	  	  	   index.index2.j,index.index2.k,index.index2.xtype,1);
			values_H1semi.insert(std::make_pair(index.index2, dd_value_x));

			return 1./std::sqrt( value_t*(value_x+dd_value_x) + pow2i<T>(2*index.index1.j)*value_t/(value_x + dd_value_x));
		}
	}
	else{
        return 1./std::sqrt( (value_x*std::pow(2.,_s*index.index2.j))
        					+ std::pow(2.,-_s*index.index1.j)*value_t*std::pow(2.,-_s*index.index2.j));
	}

}


template <typename T, typename Basis2D>
T
AdaptiveRightNormPreconditioner2D_c<T,Basis2D>::operator()(XType xtype1, int j1, long k1,
                                                   XType xtype2, int j2, long k2)
{
    return this->operator()(Index2D(Index1D(j1,k1,xtype1),Index1D(j2,k2,xtype2)));
}

}   // namespace lawa
