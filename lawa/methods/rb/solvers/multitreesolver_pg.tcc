#include <iostream>
#include <fstream>
#include <iomanip>

namespace lawa {

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
MultiTreeSolver_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, F_RHS &_F,
				TrialPrec &_trialPrec, TestPrec& _testPrec)
: trialbasis(_trialbasis), testbasis(_testbasis), Op(_Op), OpTransp(_OpTransp), F(_F),
  trialPrec(_trialPrec), testPrec(_testPrec), exact_sol(nullptr)
{}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
MultiTreeSolver_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, F_RHS &_F,
				TrialPrec &_trialPrec, TestPrec& _testPrec,
				ISWGM_Parameters& _iswgm_params, IS_Parameters& _is_params)
: iswgm_params(_iswgm_params), is_params(_is_params), trialbasis(_trialbasis), testbasis(_testbasis),
  Op(_Op), OpTransp(_OpTransp), F(_F),
  trialPrec(_trialPrec), testPrec(_testPrec), exact_sol(nullptr)
{}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
F_RHS&
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
get_rhs()
{
	return F;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
LocalOperator&
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
get_lhs()
{
	return Op;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
const TrialBasis&
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
get_trialbasis()
{
	return trialbasis;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
const TestBasis&
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
get_testbasis()
{
	return testbasis;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
TrialPrec&
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
get_trialprec()
{
	return trialPrec;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
TestPrec&
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
get_testprec()
{
	return testPrec;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
set_indexsets(const IndexSet<Index> _LambdaTrial, const IndexSet<Index> _LambdaTest)
{
	default_Lambda_trial = _LambdaTrial;
	default_Lambda_test  = _LambdaTest;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
remove_preconditioner(Coefficients<Lexicographical,T,Index> &u)
{
	for(auto& el : u){
		el.second *= trialPrec(el.first);
	}
}


template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
Coefficients<Lexicographical,typename LocalOperator::T,Index>
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
solve()
{
    IndexSet<Index> LambdaTrial = default_Lambda_trial;
    IndexSet<Index> LambdaTest = default_Lambda_test;

	Coefficients<Lexicographical,T,Index> u;
	FillWithZeros(LambdaTrial,u);

	solve(u,LambdaTrial,LambdaTest);

	return u;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
typename LocalOperator::T
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
solve(Coefficients<Lexicographical,T,Index> &u, T)
{
    IndexSet<Index> LambdaTrial, LambdaTest;
	if(u.size() > 0){

		// Take support of u as initial index set
		LambdaTrial = supp(u);

		// Check if this is the default trial index set
		bool is_default_set = true;
		if(LambdaTrial.size() != default_Lambda_trial.size()){
			is_default_set = false;
		}
		else{
			for(auto& lambda : default_Lambda_trial){
				if(LambdaTrial.find(lambda) == LambdaTrial.end()){
					is_default_set = false;
					break;
				}
			}
		}

		// If it is not the default set, we have to get some stable
		// test index set -> stable expansion
		if(!is_default_set){
			std::cerr << "MultiTreeSolver_PG: Computing stable expansion as test set. " << std::endl;
			Coefficients<Lexicographical,T,Index2D> Lambda_aux;
			getStableExpansion(trialbasis, testbasis, LambdaTrial, Lambda_aux);
			LambdaTest = supp(Lambda_aux);
		}
		// Else we just use the default test set
		else{
			LambdaTest = default_Lambda_test;
		}
	}
	else{
		LambdaTrial = default_Lambda_trial;
		LambdaTest = default_Lambda_test;
		FillWithZeros(LambdaTrial,u);
	}

	return solve(u, LambdaTrial, LambdaTest);
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
typename LocalOperator::T
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& LambdaTrial, IndexSet<Index>& LambdaTest, T)
{
    //---------------------------------------//
    //------- Initialization -----------//
    //---------------------------------------//

    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> r,s,p,Ap;
	FillWithZeros(LambdaTest,r);
	FillWithZeros(LambdaTrial,p);
	FillWithZeros(LambdaTrial,s);
	FillWithZeros(LambdaTest,Ap);

	//---------------------------------------//
	//------- CGLS  -------------------------//
	//---------------------------------------//

	// Local variables
	T alpha, beta, gamma_cgls, gamma_cgls_Prev, res_cgls;
	T dummy=0.;

	// Initial step
	Op.eval(u,r,trialPrec,testPrec);
	//r -= f;
	for(auto& lambda : LambdaTest){
		r[lambda] -= testPrec(lambda)*F(lambda);
	}
	r *= -1;
	OpTransp.eval(r,s,testPrec,trialPrec);
	p = s;
	gamma_cgls_Prev = s*s;
	iswgm_info.is_res.push_back(gamma_cgls_Prev);

	// CGLS Parameters
	T cgls_tol;
	if(is_params.adaptive_tol){
		cgls_tol = std::min(is_params.init_tol, is_params.res_reduction*gamma_cgls_Prev);
	}
	else{
		cgls_tol = is_params.absolute_tol;
	}

    if(iswgm_params.verbose){
        std::cout << "********  Starting CGLS with tolerance " << cgls_tol << " ---" << std::endl;
        std::cout << std::right;
        std::cout << "   Size of LambdaTrial " << std::setw(8) <<  LambdaTrial.size() << std::endl;
        std::cout << "   Size of LambdaTest  " << std::setw(8) << LambdaTest.size() << std::endl << std::endl;
        std::cout.precision();
        std::cout << std::left;
    }

	// CGLS Iterations
	for(std::size_t cgls_its=0; cgls_its <= is_params.max_its; ++cgls_its){

		Ap.setToZero();						// Ap = A*p
		Op.eval(p,Ap,trialPrec,testPrec);

		alpha = gamma_cgls_Prev / (Ap*Ap);
		u += alpha * p;
		r -= alpha * Ap;

		s.setToZero();
		OpTransp.eval(r,s,testPrec,trialPrec);	// s = A^T*r

		gamma_cgls = s*s;
		res_cgls = r.norm(2.);

		iswgm_info.is_resNE.push_back(gamma_cgls);
		iswgm_info.is_res.push_back(res_cgls);
		if(is_params.verbose){
			std::cout.precision(12);
			std::cout << "       CGLS Iteration " << std::setw(3) << cgls_its << ": current error NE " << std::setw(18) << sqrt(gamma_cgls)
															 << ", Au-f " << std::setw(18) << res_cgls << std::endl;
		}

		if(std::sqrt(gamma_cgls) <= cgls_tol){
			std::cerr << "       CGLS stopped with NE error " << sqrt(gamma_cgls) << ", error Au=f = " << res_cgls
				 << " after " << cgls_its << " iterations "<< std::endl;
			break;
		}

		beta  = gamma_cgls/gamma_cgls_Prev;
		p *= beta;
		p += s;
		gamma_cgls_Prev = gamma_cgls;
	}

    if(iswgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
    	if(iswgm_params.verbose){
            std::cout << "=====>  Plotting CGLS solution to file " << std::endl << std::endl;
    	}
        plot2D<T,TrialBasis,TrialPrec>(trialbasis, u, trialPrec, exact_sol, 0., 1., 0., 1., 0.01, iswgm_params.plot_filename.c_str());
    }
    
    return gamma_cgls;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
set_sol(sol_fct_2d _sol)
{
	exact_sol = _sol;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
ISWGM_Parameters&
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
access_params()
{
	return iswgm_params;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
reset_info()
{
	iswgm_info.reset();
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeSolver_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,F_RHS,TrialPrec,TestPrec>::
coarsen(Coefficients<Lexicographical,T,Index>& u){}

} // namespace lawa
