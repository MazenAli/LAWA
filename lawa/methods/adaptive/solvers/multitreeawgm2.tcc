#include <lawa/methods/adaptive/algorithms/algorithms.h>

#include <ios>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace lawa {

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM2(const Basis &_basis, LocalOperator &_Op, RHS &_F, Preconditioner &_Prec)
: basis(_basis), Op(_Op), F(_F), Prec(_Prec), exact_sol(nullptr)
{}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM2(const Basis &_basis, LocalOperator &_Op, RHS &_F, Preconditioner &_Prec,
				AWGM_Parameters& _awgm_params, IS_Parameters& _is_params)
: awgm_params(_awgm_params), is_params(_is_params),
  basis(_basis), Op(_Op), F(_F), Prec(_Prec), exact_sol(nullptr)
{}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
set_initial_indexset(const IndexSet<Index> _LambdaTrial)
{
	default_init_Lambda = _LambdaTrial;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
remove_preconditioner(Coefficients<Lexicographical,T,Index> &u)
{
	for(auto& el : u){
		el.second *= Prec(el.first);
	}
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
RHS&
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
get_rhs()
{
	return F;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
LocalOperator&
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
get_lhs()
{
	return Op;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
const Basis&
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
get_trialbasis()
{
	return basis;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
const Basis&
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
get_testbasis()
{
	return basis;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
Preconditioner&
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
get_trialprec()
{
	return Prec;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
Preconditioner&
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
get_testprec()
{
	return Prec;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
Coefficients<Lexicographical,typename LocalOperator::T,Index>
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
solve()
{
    IndexSet<Index> Lambda = default_init_Lambda;

	Coefficients<Lexicographical,T,Index> u;
	FillWithZeros(Lambda,u);

	solve(u,Lambda);

	return u;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
typename LocalOperator::T
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
solve(Coefficients<Lexicographical,T,Index> &u, T old_res_norm)
{
    IndexSet<Index> Lambda;
	if(u.size() > 0){
		// Take support of u as initial index set
		Lambda = supp(u);
	}
	// Else we just use the default test set
	else{
		Lambda = default_init_Lambda;
		FillWithZeros(Lambda,u);
	}

	return solve(u, Lambda, old_res_norm);

}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
typename LocalOperator::T
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& Lambda, T old_res_norm)
{
    //---------------------------------------//
    //------- AWGM Initialization -----------//
    //---------------------------------------//

    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> r(awgm_params.hashmapsize),
                                            p(awgm_params.hashmapsize),
                                            Ap(awgm_params.hashmapsize),
                                            res(awgm_params.hashmapsize),     // approximate residual for f-Au
                                            u_leafs(awgm_params.hashmapsize); // "leafs" of u

	FillWithZeros(Lambda,u_leafs);
	FillWithZeros(Lambda,res);

    // Default for initial cg tolerance computation if adaptive
    T res_norm = old_res_norm;

    //---------------------------------------//
    //------- AWGM Iterations ---------------//
    //---------------------------------------//
    for(std::size_t awgm_its = 0; awgm_its <= awgm_params.max_its; ++awgm_its){

    	if (Lambda.size()>awgm_params.max_basissize){
    	    std::cerr << "AWGM reached maximal basis size " << awgm_params.max_basissize << ": "
    	    		  << "Residual = " << res_norm << " "
    	              << ", awgm_tol = " << awgm_params.tol << std::endl << std::endl;
    	    break;
    	}

        awgm_info.sizeLambda.push_back(Lambda.size());

        //---------------------------------------//
        //------- CGLS  -------------------------//
        //---------------------------------------//

		// Re-initialize vectors from last iteration
        r.clear();
        p.clear();
        Ap.clear();
		FillWithZeros(Lambda,r);
		FillWithZeros(Lambda,p);
		FillWithZeros(Lambda,Ap);

		// CG Parameters
		T cg_tol;
		if(is_params.adaptive_tol){
			cg_tol = std::min(is_params.init_tol, is_params.res_reduction*res_norm);
		}
		else{
			cg_tol = is_params.absolute_tol;
		}

        if(awgm_params.verbose){
            std::cout << "******** Iteration " << awgm_its << std::endl;
            std::cout << std::right;
            std::cout << "   Current size of Lambda " << std::setw(8) <<  Lambda.size() << std::endl;
            std::cout.precision();
            std::cout << std::left;
            std::cout << "   --- Starting CG with tolerance " << cg_tol << " ---" << std::endl;
        }

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

		// CG Iterations
		for(std::size_t cg_its=0; cg_its <= is_params.max_its; ++cg_its){

			Ap.setToZero();						// Ap = A*p
			Op.eval(p,Ap,Prec);

			alpha = res_cg_prev / (p*Ap);
			u += alpha * p;
			r -= alpha * Ap;

			res_cg = r*r;

			if(is_params.verbose){
				std::cout.precision(12);
	            std::cout << "       CG Iteration " << std::setw(3) << cg_its << ": current error " << std::setw(18) << sqrt(res_cg) << std::endl;
			}

			if(std::sqrt(res_cg) <= cg_tol){
                std::cerr << "       CG stopped with error " << sqrt(res_cg) << " after "
                		  << cg_its << " iterations "<< std::endl;
                awgm_info.cg_its.push_back(cg_its);
				break;
			}

			beta  = res_cg/res_cg_prev;
			p *= beta;
			p += r;
			res_cg_prev = res_cg;
		}

		if(awgm_info.cg_its.size() < awgm_info.sizeLambda.size()){
            awgm_info.cg_its.push_back(is_params.max_its);
		}

        //---------------------------------------//
        //---- COMPUTE APPROXIMATE RESIDUAL -----//
        //---------------------------------------//

		// Compute cone around solution u
		IndexSet<Index2D>	Cdiff_u_leafs;	// new indizes in cone
        res.setToZero();
        extendMultiTree(basis, u_leafs, res, Cdiff_u_leafs, "standard", false, true);

        // Compute rhs on expanded test index set
        IndexSet<Index2D> LambdaRes = supp(res);

		// Compute residual Au - f on expanded test index set
        res.setToZero();
        Op.eval(u,res,Prec); 	// res = A*u
        //res -= f;
//	    for(auto& lambda : LambdaRes){
//	    	res[lambda] -= Prec(lambda)*F(lambda);
//	    }
		{
			Coefficients<Lexicographical,T,Index> f = F(LambdaRes);
		    for(auto& lambda : LambdaRes){
		    	res[lambda] -= Prec(lambda)*f[lambda];
		    }
		}

        res_norm = res.norm(2.);

        awgm_info.awgm_res.push_back(res_norm);
        awgm_info.sizeLambdaRes.push_back(LambdaRes.size());

        if(awgm_params.verbose){
            std::cout << "   --- Approximate Residual ---" << std::endl;
            std::cout << "       U_leafs:        " << std::setw(20) << u_leafs.size() << std::endl;
            std::cout << "       Residual:  " << std::setw(20) << res_norm
            								  << " (on " << std::setw(8) << std::right <<  LambdaRes.size()
            								  << std::left << " indizes)" << std::endl;
        }

        // Test for convergence
        if(res_norm <= awgm_params.tol){
            std::cerr << "AWGM target tolerance reached after " << awgm_its << " iterations: "
            		  << "Residual = " << res_norm << " "
                      << ", awgm_tol = " << awgm_params.tol << std::endl << std::endl;

            if(awgm_params.print_info){
            	if(awgm_params.verbose){
                    std::cout << "=====>  Writing AWGM CG information to file " << std::endl << std::endl;
            	}
            	awgm_info.print(awgm_params.info_filename.c_str());
            }

            if(awgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
            	if(awgm_params.verbose){
                    std::cout << "=====>  Plotting AWGM CG solution to file " << std::endl << std::endl;
            	}
                plot2D<T,Basis,Preconditioner>(basis, u, Prec, exact_sol, 0., 1., 0., 1., 0.01, awgm_params.plot_filename.c_str());
            }

            if(awgm_params.clear_solver){
                Op.clear();
                F.clear();
            }
            
            return res_norm;
        }

        //---------------------------------------//
        //---- COMPUTE NEXT INDEX SET -----------//
        //---------------------------------------//

        // Remove indizes from last iteration
        T P_Lambda_Res_square = 0.;
		for(auto& lambda : Lambda){
			P_Lambda_Res_square += std::pow(r[lambda], (T)2.);
			res.erase(lambda);
		}

		// Compute buckets
		T threshbound = std::sqrt(1.-awgm_params.alpha*awgm_params.alpha) * res.norm((T)2.)/std::sqrt(T(res.size()));
		Coefficients<Bucket,T,Index2D> r_bucket;
		r_bucket.bucketsort(res, threshbound);

        if(awgm_params.verbose){
            std::cout << "   --- Computing next index set ---" << std::endl;
            std::cout << "       Threshbound:         " << threshbound <<  std::endl;
            std::cout << "       ||P_{Lambda}r ||_2:  " << std::sqrt(P_Lambda_Res_square) << std::endl;
            std::cout << "       alpha*Residual:   " << awgm_params.alpha*res_norm << std::endl;
        }

		// Add buckets to dummy vector: p
		for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
			P_Lambda_Res_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
			r_bucket.addBucketToCoefficients(p,i);

			if(awgm_params.verbose){
				std::cerr << "          Bucket " << i << ": L2-norm " << r_bucket.bucket_ell2norms[i] << std::endl;
			}

			if (P_Lambda_Res_square >= awgm_params.alpha*res_norm*awgm_params.alpha*res_norm) {
				break;
			}
		}

		// Replace indizes from last iteration
		for(auto& lambda : Lambda){
			res.insert(std::make_pair<decltype(lambda),T>(lambda,0.));
		}

		// New LambdaTrial: Add new indizes to u, complete to multitree
        IndexSet<Index2D> LambdaDiff;
        for(auto& coeff : p){
			if(u.find(coeff.first) == u.end()){
				completeMultiTree(basis, coeff.first, u, LambdaDiff, 0, true);
			}
		}

        // Compute new u_leafs
        u_leafs.clear();
        FillWithZeros(LambdaDiff,u_leafs);

        Lambda = supp(u);

        if(awgm_params.verbose){
            std::cout << "       Lambda:  raw extension  " <<  std::setw(8) << std::right <<  p.size() << std::endl;
            std::cout << "                multitree ext. " <<  std::setw(8) <<  LambdaDiff.size() << std::endl;
            std::cout << "                total size     " <<  std::setw(8) <<  Lambda.size() << std::endl;
            if(awgm_params.verbose_extra){
            	std::cout << LambdaDiff;
            }
            std::cout << std::endl;
        }
    }

    std::cerr << "AWGM reached maximal iteration number " << awgm_params.max_its << ": "
    		  << "Residual = " << res_norm << " "
              << ", awgm_tol = " << awgm_params.tol << std::endl << std::endl;

    if(awgm_params.print_info){
    	if(awgm_params.verbose){
            std::cout << "=====>  Writing AWGM CG information to file " << std::endl << std::endl;
    	}
    	awgm_info.print(awgm_params.info_filename.c_str());
    }
    if(awgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
    	if(awgm_params.verbose){
            std::cout << "=====>  Plotting AWGM CG solution to file " << std::endl << std::endl;
    	}
        plot2D<T,Basis,Preconditioner>(basis, u, Prec, exact_sol, 0., 1., 0., 1., 0.01, awgm_params.plot_filename.c_str());
    }

    if(awgm_params.clear_solver){
        Op.clear();
        F.clear();
    }
    
    return res_norm;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
coarsen(Coefficients<Lexicographical,T,Index>& u, T delta)
{
	if(u.size() > 0){
		// Take support of u as initial index set
	    IndexSet<Index> Lambda = supp(u);

	    Coefficients<Lexicographical,T,Index2D> res(awgm_params.hashmapsize);
		FillWithZeros(Lambda,res);

		// res = Au - f
		Op.eval(u,res,Prec);
		{
			Coefficients<Lexicographical,T,Index> f = F(Lambda);
		    for(auto& lambda : Lambda){
		    	res[lambda] -= Prec(lambda)*f[lambda];
		    }
		}
        T res_norm = res.norm(2.);

		// Compute buckets
		T threshbound = std::sqrt(1.-delta*delta) * res_norm/std::sqrt(T(res.size()));
		Coefficients<Bucket,T,Index2D> r_bucket;
		r_bucket.bucketsort(res, threshbound);

        if(awgm_params.verbose){
            std::cout << "   --- Coarsen solution ---" << std::endl;
            std::cout << "       Threshbound:         " << threshbound <<  std::endl;
            std::cout << "       alpha*Residual:   " << awgm_params.alpha*res_norm << std::endl;
        }

		// Reconstruct u by adding buckets
		u.clear();
		T P_Lambda_Res_square = 0;
		for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
			P_Lambda_Res_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
			r_bucket.addBucketToCoefficients(u,i);

			if(awgm_params.verbose){
				std::cerr << "          Bucket " << i << ": L2-norm " << r_bucket.bucket_ell2norms[i] << std::endl;
			}

			if (P_Lambda_Res_square >= delta*res_norm*delta*res_norm) {
				break;
			}
		}
	}
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
coarsen(Coefficients<Lexicographical,T,Index>& u)
{
	coarsen(u, awgm_params.alpha);
}



template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
typename LocalOperator::T
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
estimate_residual_norm(Coefficients<Lexicographical,T,Index>& u)
{
    //---------------------------------------//
    //---- COMPUTE APPROXIMATE RESIDUAL -----//
    //---------------------------------------//
    Coefficients<Lexicographical,T,Index2D> u_copy(u);
    
    // Precondition u 
    for(auto& el : u_copy){
		el.second /= Prec(el.first);
	}

	// Compute cone around solution u
	Coefficients<Lexicographical,T,Index2D> res(awgm_params.hashmapsize);
    FillWithZeros(supp(u_copy),res);
    extendMultiTree(basis, u_copy, res, "standard", false, true);

    // Compute rhs on expanded test index set
    IndexSet<Index2D> LambdaRes = supp(res); 
    std::cout << "Estimate residual norm on " << LambdaRes.size() << " indizes." << std::endl;

	// Compute residual Au - f on expanded test index set
    res.setToZero();
    Op.eval(u_copy,res,Prec); 	// res = A*u
    //res -= f;
	Coefficients<Lexicographical,T,Index> f = F(LambdaRes);
	for(auto& lambda : LambdaRes){
		res[lambda] -= Prec(lambda)*f[lambda];
	}

    return res.norm(2.);
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
set_sol(sol_fct_2d _sol)
{
	exact_sol = _sol;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
AWGM_Parameters&
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
access_params()
{
	return awgm_params;
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
reset_info()
{
	return awgm_info.reset();
}

} // namespace lawa
