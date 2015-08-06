#include <cmath>

namespace lawa {

template <typename T, typename Basis2D>
AdaptiveLeftNormPreconditioner2D<T,Basis2D>::AdaptiveLeftNormPreconditioner2D(const Basis2D &basis, T s)
    : _basis(basis), _s(s),
      _integral(basis.second, basis.second)
{
}

/* Do _all_ calculations here (and only redirect from the operator()(j1,k1,xtype1,..):
 * 		Usually, we call the preconditioner in an adaptive setting using the Index2D, so
 * 		that this saves one function call.
 */
template <typename T, typename Basis2D>
T
AdaptiveLeftNormPreconditioner2D<T,Basis2D>::operator()(const Index2D &index)
{
	// See if we have calculated the entry before
	auto it = values_L2.find(index.index2);

	// If yes, reuse the values
	if( it != values_L2.end()){
		if(_s==2.){
			auto it_H1 = values_H1semi.find(index.index2);
			assert(it_H1 != values_H1semi.end());
			return 1./std::sqrt((*it).second + (*it_H1).second);
		}
		else{
			return 1./std::sqrt((*it).second + std::pow(2.,_s*index.index2.j));
		}
	}
	// Else compute and store
	else{

	    T value = _integral(index.index2.j,index.index2.k,index.index2.xtype,0,index.index2.j,index.index2.k,index.index2.xtype,0);
	    values_L2.insert(std::make_pair(index.index2,value));

	    if (_s==2.) {
	        // Calculate H1-Norm of Basis Function using Integrals
	    	T value_deriv = _integral(index.index2.j,index.index2.k,index.index2.xtype,1,index.index2.j,index.index2.k,index.index2.xtype,1);
	        values_H1semi.insert(std::make_pair(index.index2,value_deriv));
	        return 1./std::sqrt(value + value_deriv);
	    } else {
	        // Calculate H1-Norm of Basis Function using Scaling of the L2-norm (assumed to be equivalent to 1)
	        return 1./std::sqrt(value + std::pow(2.,_s*index.index2.j));
	    }
	}
}

template <typename T, typename Basis2D>
T
AdaptiveLeftNormPreconditioner2D<T,Basis2D>::operator()(XType xtype1, int j1, long k1,
                                                		XType xtype2, int j2, long k2)
{
    return this->operator()(Index2D(Index1D(j1,k1,xtype1),Index1D(j2,k2,xtype2)));
}

}   // namespace lawa

