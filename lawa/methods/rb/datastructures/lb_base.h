/*
 * lb_base.h
 *
 *  Created on: 09.04.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_LB_BASE_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_LB_BASE_H_

#include <vector>

#include <extensions/flens/flens.h>
#include <lawa/methods/rb/datastructures/thetastructure.h>


namespace lawa {

/**
 * Class for calculations of coercivity & inf-sup constants etc.
 */

template <typename ParamType, typename Truth>
class LB_Base {

	typedef typename Truth::T 													T;
    typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >    		SparseMatrixT;
    typedef flens::SparseSyMatrix<flens::extensions::CRS<T,flens::CRS_UpperTriangular> >    SparseSymMatrixT;
    typedef flens::DenseVector<flens::Array<T> >                        		DenseVectorT;
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  		FullColMatrixT;

public:

    LB_Base(Truth& _truth, ThetaStructure<ParamType>& _thetas_lhs, bool _preconditioned_truth = false);

    T
    calculate_alpha(ParamType& mu);

    template <typename IndexSetType>
    void
    assemble_matrices_for_alpha_computation(IndexSetType& indexset);

private:

    Truth& 							truth;
    ThetaStructure<ParamType>& 		thetas;

    bool 							preconditioned_truth;

	std::vector<SparseMatrixT> 		A_u_u_matrices;
	SparseMatrixT					I_Y_u_u_matrix;

};


} // namespace lawa

#include <lawa/methods/rb/datastructures/lb_base.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_LB_BASE_H_ */
