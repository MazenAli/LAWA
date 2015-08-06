#include <lawa/methods/adaptive/algorithms/algorithms.h>

#include <ios>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace lawa {

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
MultiTreeSolver(const Basis &_basis, LocalOperator &_Op, F_RHS &_F, Preconditioner &_Prec)
: basis(_basis), Op(_Op), F(_F), Prec(_Prec), exact_sol(nullptr)
{}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
MultiTreeSolver(const Basis &_basis, LocalOperator &_Op, F_RHS &_F, Preconditioner &_Prec,
				ISWGM_Parameters& _iswgm_params, IS_Parameters& _is_params)
: iswgm_params(_iswgm_params), is_params(_is_params),
  basis(_basis), Op(_Op), F(_F), Prec(_Prec), exact_sol(nullptr)
{}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
void
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
set_indexset(const IndexSet<Index> _LambdaTrial)
{
	default_Lambda = _LambdaTrial;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
void
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
remove_preconditioner(Coefficients<Lexicographical,T,Index> &u)
{
	for(auto& el : u){
		el.second *= Prec(el.first);
	}
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
F_RHS&
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
get_rhs()
{
	return F;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
LocalOperator&
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
get_lhs()
{
	return Op;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
const Basis&
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
get_trialbasis()
{
	return basis;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
const Basis&
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
get_testbasis()
{
	return basis;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
Preconditioner&
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
get_trialprec()
{
	return Prec;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
Preconditioner&
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
get_testprec()
{
	return Prec;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
Coefficients<Lexicographical,typename LocalOperator::T,Index>
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
solve()
{
    IndexSet<Index> Lambda = default_Lambda;

	Coefficients<Lexicographical,T,Index> u;
	FillWithZeros(Lambda,u);

	solve(u,Lambda);

	return u;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
typename LocalOperator::T
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
solve(Coefficients<Lexicographical,T,Index> &u, T)
{
    IndexSet<Index> Lambda;
	if(u.size() > 0){
		// Take support of u as initial index set
		Lambda = supp(u);
	}
	// Else we just use the default test set
	else{
		Lambda = default_Lambda;
		FillWithZeros(Lambda,u);
	}

	return solve(u, Lambda);

}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
typename LocalOperator::T
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& Lambda, T)
{
    //---------------------------------------//
    //------- Initialization -----------//
    //---------------------------------------//

    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> r, p, Ap;

	FillWithZeros(Lambda,r);
	FillWithZeros(Lambda,p);
	FillWithZeros(Lambda,Ap);

	//---------------------------------------//
	//------- CG  ---------------------------//
	//---------------------------------------//

	// Local variables
	T alpha, beta, res_cg, res_cg_prev;
	T dummy=0.;

	// Initial step
	Op.eval(u,r,Prec);
	//r -= f;
	{
		Coefficients<Lexicographical,T,Index> f = F(Lambda);
		for(auto& lambda : Lambda){
			r[lambda] -= Prec(lambda)*f[lambda];
		}
	}

	r *= -1;
	p = r;
	res_cg_prev = r*r;

	// CG Parameters
	T cg_tol;
	if(is_params.adaptive_tol){
		cg_tol = std::min(is_params.init_tol, is_params.res_reduction*res_cg_prev);
	}
	else{
		cg_tol = is_params.absolute_tol;
	}

    if(iswgm_params.verbose){
        std::cout << "********  Starting CG with tolerance " << cg_tol << " ---" << std::endl;
        std::cout << std::right;
        std::cout << "   Size of Lambda " << std::setw(8) <<  Lambda.size() << std::endl;
        std::cout.precision();
        std::cout << std::left;
    }
	iswgm_info.is_res.push_back(res_cg_prev);


	// CG Iterations
	for(std::size_t cg_its=0; cg_its <= is_params.max_its; ++cg_its){

		Ap.setToZero();						// Ap = A*p
		Op.eval(p,Ap,Prec);

		alpha = res_cg_prev / (p*Ap);
		u += alpha * p;
		r -= alpha * Ap;

		res_cg = r*r;

		iswgm_info.is_res.push_back(res_cg);
		if(is_params.verbose){
			std::cout.precision(12);
			std::cout << "       CG Iteration " << std::setw(3) << cg_its << ": current error " << std::setw(18) << sqrt(res_cg) << std::endl;
		}

		if(std::sqrt(res_cg) <= cg_tol){
			std::cerr << "       CG stopped with error " << sqrt(res_cg) << " after "
					  << cg_its << " iterations "<< std::endl;
			break;
		}

		beta  = res_cg/res_cg_prev;
		p *= beta;
		p += r;
		res_cg_prev = res_cg;
	}

	if(iswgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
		if(iswgm_params.verbose){
			std::cout << "=====>  Plotting CG solution to file " << std::endl << std::endl;
		}
		plot2D<T,Basis,Preconditioner>(basis, u, Prec, exact_sol, 0., 1., 0., 1., 0.01, iswgm_params.plot_filename.c_str());
	}

    return res_cg;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
void
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
set_sol(sol_fct_2d _sol)
{
	exact_sol = _sol;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
ISWGM_Parameters&
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
access_params()
{
	return iswgm_params;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
void
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
reset_info()
{
	iswgm_info.reset();
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
void
MultiTreeSolver<Index,Basis,LocalOperator,F_RHS,Preconditioner>::
coarsen(Coefficients<Lexicographical,T,Index>& u){}


} // namespace lawa
