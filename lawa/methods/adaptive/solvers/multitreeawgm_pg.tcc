#include <lawa/methods/adaptive/algorithms/algorithms.h>

#include <ios>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace lawa {

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
MultiTreeAWGM_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, RHS &_F,
				TrialPrec &_trialPrec, TestPrec& _testPrec)
: trialbasis(_trialbasis), testbasis(_testbasis), Op(_Op), OpTransp(_OpTransp), F(_F),
  trialPrec(_trialPrec), testPrec(_testPrec), exact_sol(nullptr)
{}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
MultiTreeAWGM_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, RHS &_F,
				TrialPrec &_trialPrec, TestPrec& _testPrec,
				AWGM_PG_Parameters& _awgm_params, IS_Parameters& _is_params)
: awgm_params(_awgm_params), is_params(_is_params),
  trialbasis(_trialbasis), testbasis(_testbasis), Op(_Op), OpTransp(_OpTransp), F(_F),
  trialPrec(_trialPrec), testPrec(_testPrec), exact_sol(nullptr)
{}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
RHS&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_rhs()
{
	return F;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
LocalOperator&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_lhs()
{
	return Op;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
const TrialBasis&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_trialbasis()
{
	return trialbasis;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
const TestBasis&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_testbasis()
{
	return testbasis;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
TrialPrec&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_trialprec()
{
	return trialPrec;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
TestPrec&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_testprec()
{
	return testPrec;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
set_initial_indexsets(const IndexSet<Index> _LambdaTrial, const IndexSet<Index> _LambdaTest)
{
	default_init_Lambda_trial = _LambdaTrial;
	default_init_Lambda_test  = _LambdaTest;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
remove_preconditioner(Coefficients<Lexicographical,T,Index> &u)
{
	for(auto& el : u){
		el.second *= trialPrec(el.first);
	}
}


template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
Coefficients<Lexicographical,typename LocalOperator::T,Index>
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
solve()
{
    IndexSet<Index> LambdaTrial = default_init_Lambda_trial;
    IndexSet<Index> LambdaTest = default_init_Lambda_test;

	Coefficients<Lexicographical,T,Index> u;
	FillWithZeros(LambdaTrial,u);

	solve(u,LambdaTrial,LambdaTest);

	return u;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
solve(Coefficients<Lexicographical,T,Index> &u)
{
    IndexSet<Index> LambdaTrial, LambdaTest;
	if(u.size() > 0){

		// Take support of u as initial index set
		LambdaTrial = get_indexset(u);

		// Check if this is the default trial index set
		bool is_default_set = true;
		if(LambdaTrial.size() != default_init_Lambda_trial.size()){
			is_default_set = false;
		}
		else{
			for(auto& lambda : default_init_Lambda_trial){
				if(LambdaTrial.find(lambda) == LambdaTrial.end()){
					is_default_set = false;
					break;
				}
			}
		}

		// If it is not the default set, we have to get some stable
		// test index set -> stable expansion
		if(!is_default_set){
			std::cerr << "Computing stable expansion as test set. " << std::endl;
			Coefficients<Lexicographical,T,Index2D> Lambda_aux;
	        switch(awgm_params.stable_exp_u){
				case FullExpansion:
					getStableExpansion(trialbasis, testbasis, LambdaTrial, Lambda_aux);
					break;
				case WoMixedHWExpansion:
					getStableExpansion_woMixedHW(trialbasis, testbasis, LambdaTrial, Lambda_aux);
					break;
				case OnlyTemporalHWExpansion:
					getStableExpansion_onlyTemporalHW(trialbasis, testbasis, LambdaTrial, Lambda_aux);
					break;
				default:
					std::cerr << "Stable Expansion type doesn't exist!" << std::endl;
					break;
	        }
			LambdaTest = get_indexset(Lambda_aux);
		}
		// Else we just use the default test set
		else{
			LambdaTest = default_init_Lambda_test;
		}
	}
	else{
		LambdaTrial = default_init_Lambda_trial;
		LambdaTest = default_init_Lambda_test;
		FillWithZeros(LambdaTrial,u);
	}

	solve(u, LambdaTrial, LambdaTest);
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& LambdaTrial, IndexSet<Index>& LambdaTest)
{
    //---------------------------------------//
    //------- AWGM Initialization -----------//
    //---------------------------------------//

    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> r(awgm_params.hashmapsize_test),
    										s(awgm_params.hashmapsize_trial),
                                            p(awgm_params.hashmapsize_trial),
                                            Ap(awgm_params.hashmapsize_test),
                                            res(awgm_params.hashmapsize_test),     // approximate residual for f-Au
                                    		resNE(awgm_params.hashmapsize_trial),   // cone around u_leafs
                                            u_leafs(awgm_params.hashmapsize_trial), // "leafs" of u
                                            resNE_tmp(awgm_params.hashmapsize_trial); // cone around u_leafs if we further extend resNE


    if(u.size() != LambdaTrial.size()){
    	FillWithZeros(LambdaTrial, u);
    }
	FillWithZeros(LambdaTrial,u_leafs);
	FillWithZeros(LambdaTest,res);
	FillWithZeros(LambdaTest,Ap); // our dummy vector for the test index set extension

    // Default for initial cgls tolerance computation if adaptive
    T resNE_norm = 1.;
    T res_norm = 1.;

    //---------------------------------------//
    //------- AWGM Iterations ---------------//
    //---------------------------------------//
    for(std::size_t awgm_its = 0; awgm_its <= awgm_params.max_its; ++awgm_its){

    	if (u.size()>awgm_params.max_basissize){
    		break;
    	}

        awgm_info.sizeLambdaTrial.push_back(u.size());
        awgm_info.sizeLambdaTest.push_back(Ap.size());

        //---------------------------------------//
        //------- CGLS  -------------------------//
        //---------------------------------------//

		// Re-initialize vectors from last iteration
        // We assume that the indexsets never shrink
        for(auto& entry : u){
        	p[entry.first] = 0.;
        	s[entry.first] = 0.;
        }
        for(auto& entry: Ap){
        	r[entry.first] = 0.;
        	entry.second = 0.;
        }

		// CGLS Parameters
		T cgls_tol;
		if(is_params.adaptive_tol){
			cgls_tol = std::min(is_params.init_tol, is_params.res_reduction*resNE_norm);
		}
		else{
			cgls_tol = is_params.absolute_tol;
		}

        if(awgm_params.verbose){
            std::cout << "******** Iteration " << awgm_its << std::endl;
            std::cout << std::right;
            std::cout << "   Current size of LambdaTrial " << std::setw(8) <<  u.size() << std::endl;
            std::cout << "   Current size of LambdaTest  " << std::setw(8) << Ap.size() << std::endl << std::endl;
            std::cout.precision();
            std::cout << std::left;
            std::cout << "   --- Starting CGLS with tolerance " << cgls_tol << " ---" << std::endl;
        }

		// Local variables
		T alpha, beta, gamma_cgls, gamma_cgls_Prev, res_cgls;
		T dummy=0.;

		// Initial step
		Op.eval(u,r,trialPrec,testPrec);
		//r -= f;
	    for(auto& entry : r){
	    	entry.second -= testPrec(entry.first)*F(entry.first);
	    }
		r *= -1;
		OpTransp.eval(r,s,testPrec,trialPrec);
		p = s;
		gamma_cgls_Prev = s*s;

		// CGLS Iterations
		for(std::size_t cgls_its=0; cgls_its <= is_params.max_its; ++cgls_its){

			Ap.setToZero();// Ap = A*p
			//std::cout << "Eval Ap = A * p" << std::endl;
			Op.eval(p,Ap,trialPrec,testPrec);

			alpha = gamma_cgls_Prev / (Ap*Ap);
			u += alpha * p;
			r -= alpha * Ap;

			s.setToZero();

			//std::cout << "Eval s = A^T * r" << std::endl;
			OpTransp.eval(r,s,testPrec,trialPrec);	// s = A^T*r

			gamma_cgls = s*s;
			res_cgls = r.norm(2.);

			if(is_params.verbose){
				std::cout.precision(12);
	            std::cout << "       CGLS Iteration " << std::setw(3) << cgls_its << ": current error NE " << std::setw(18) << sqrt(gamma_cgls)
	            												 << ", Au-f " << std::setw(18) << res_cgls << std::endl;
			}

			if(std::sqrt(gamma_cgls) <= cgls_tol){
                std::cerr << "       CGLS stopped with NE error " << sqrt(gamma_cgls) << ", error Au=f = " << res_cgls
                	 << " after " << cgls_its << " iterations "<< std::endl;
                awgm_info.cgls_its.push_back(cgls_its);
				break;
			}

			beta  = gamma_cgls/gamma_cgls_Prev;
			p *= beta;
			p += s;
			gamma_cgls_Prev = gamma_cgls;
		}

		if(awgm_info.cgls_its.size() < awgm_info.sizeLambdaTrial.size()){
            awgm_info.cgls_its.push_back(is_params.max_its);
		}

    	if(awgm_params.write_intermediary_solutions){
    		std::stringstream filename;
    		filename << awgm_params.intermediary_solutions_filename << "_" << awgm_its << ".txt";
    	    saveCoeffVector2D(u, trialbasis, filename.str().c_str());
    	}


        //---------------------------------------//
        //---- COMPUTE APPROXIMATE RESIDUAL -----//
        //---------------------------------------//

		// Compute cone around solution u
		IndexSet<Index2D>	Cdiff_u_leafs;	// new indizes in cone

        if(awgm_params.reset_res){
        	res.clear();
        }
        else{
            res.setToZero();
        }

        // If we have more than one cone extension for resNE, use auxiliary index set resNE_tmp
        // for the cone around u_leafs, as it blows up otherwise (we never "reset" resNE)
        if(awgm_params.res_construction == SimpleStableExpansion || awgm_params.res_construction == ParallelMTCones){
            extendMultiTree(trialbasis, u_leafs, resNE, Cdiff_u_leafs, "standard", false, true);

            // Compute stable expansion in test basis to get res
            if(awgm_params.res_construction != ParallelMTCones){
    			switch(awgm_params.stable_exp_res){
    				case FullExpansion:
    					getStableExpansion(trialbasis, testbasis, resNE, res);
    					break;
    				case WoMixedHWExpansion:
    					getStableExpansion_woMixedHW(trialbasis, testbasis, resNE, res);
    					break;
    				case OnlyTemporalHWExpansion:
    					getStableExpansion_onlyTemporalHW(trialbasis, testbasis, resNE, res);
    					break;
    				default:
    					std::cerr << "Stable Expansion type doesn't exist!" << std::endl;
    					break;
    			}
            }
            else{
            	extendMultiTree(testbasis, Ap, res, "standard", false, true);
            }
        }
        else{
            extendMultiTree(trialbasis, u_leafs, resNE_tmp, Cdiff_u_leafs, "standard", false, true);

            // Compute stable expansion of resNE_tmp in test basis to get res
			switch(awgm_params.stable_exp_res){
				case FullExpansion:
					getStableExpansion(trialbasis, testbasis, resNE_tmp, res);
					break;
				case WoMixedHWExpansion:
					getStableExpansion_woMixedHW(trialbasis, testbasis, resNE_tmp, res);
					break;
				case OnlyTemporalHWExpansion:
					getStableExpansion_onlyTemporalHW(trialbasis, testbasis, resNE_tmp, res);
					break;
				default:
					std::cerr << "Stable Expansion type doesn't exist!" << std::endl;
					break;
			}

			if(awgm_params.res_construction == DoubleStableExpansion){
				resNE.clear();
				getStableExpansion(testbasis, trialbasis, res, resNE);
			}
			if(awgm_params.res_construction == DoubleMTExtension){
		        extendMultiTree(trialbasis, resNE_tmp, resNE, Cdiff_u_leafs, "standard", false, true);
			}
        }

        // Compute rhs on expanded test index set

		// Compute residual Au - f on expanded test index set
        res.setToZero();
        Op.eval(u,res,trialPrec,testPrec); 	// res = A*u
        //res -= f;
		for(auto& entry : res){
			entry.second -= testPrec(entry.first)*F(entry.first);
		}
        res_norm = res.norm(2.);

        awgm_info.awgm_res.push_back(res_norm);
        awgm_info.sizeLambdaRes.push_back(res.size());

        // Compute residual of NE: A^TA*u - A^T*f
        resNE.setToZero();
        OpTransp.eval(res,resNE,testPrec,trialPrec); 	// resNE = A^T*res
        resNE_norm = resNE.norm(2.);

        awgm_info.awgm_resNE.push_back(resNE_norm);
        awgm_info.sizeLambdaResNE.push_back(resNE.size());

        if(awgm_params.verbose){
            std::cout << "   --- Approximate Residual ---" << std::endl;
            std::cout << "       U_leafs:        " << std::setw(20) << u_leafs.size() << std::endl;
            std::cout << "       Residual Au-f:  " << std::setw(20) << res_norm
            									   << " (on " << std::setw(8) << std::right <<  res.size() << std::left << " indizes)" << std::endl;
            std::cout << "       Residual NE  :  " << std::setw(20) << resNE_norm
            									   << " (on " << std::setw(8) << std::right << resNE.size()<< std::left << " indizes)" << std::endl;
        }

        if(awgm_params.write_intermediary_solutions && awgm_params.print_info){
        	awgm_info.print(awgm_params.info_filename.c_str());
        }

        // Test for convergence
        if(res_norm <= awgm_params.tol_primal || resNE_norm <= awgm_params.tol){
            std::cerr << "AWGM target tolerance reached after " << awgm_its << " iterations: "
            		  << "Residual NE = " << resNE_norm << " "
            		  << ", Residual Au-f = " << res_norm << " "
                      << ", awgm_tol = " << awgm_params.tol << " "
                      << ", primal_tol = " << awgm_params.tol_primal << std::endl<< std::endl;

            if(awgm_params.print_info){
            	if(awgm_params.verbose){
                    std::cout << "=====>  Writing AWGM CGLS information to file " << std::endl << std::endl;
            	}
            	awgm_info.print(awgm_params.info_filename.c_str());
            }

            if(awgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
            	if(awgm_params.verbose){
                    std::cout << "=====>  Plotting AWGM CGLS solution to file " << std::endl << std::endl;
            	}
                plot2D<T,TrialBasis,TrialPrec>(trialbasis, u, trialPrec, exact_sol, 0., 1., 0., 1., 0.01, awgm_params.plot_filename.c_str());
            }
            
            if(awgm_params.clear_solver){
                Op.clear();
                OpTransp.clear();
                F.clear();
            }

            return;
        }

        //---------------------------------------//
        //---- COMPUTE NEXT INDEX SET -----------//
        //---------------------------------------//

        // Remove indizes from last iteration
        T P_Lambda_ResNE_square = 0.;
		for(auto& entry : s){
			P_Lambda_ResNE_square += std::pow(entry.second, (T)2.);
			resNE.erase(entry.first);
		}

		// Compute buckets
		T threshbound = std::sqrt(1.-awgm_params.alpha*awgm_params.alpha) * resNE.norm((T)2.)/std::sqrt(T(resNE.size()));
		Coefficients<Bucket,T,Index2D> r_bucket;
		r_bucket.bucketsort(resNE, threshbound);

        if(awgm_params.verbose){
            std::cout << "   --- Computing next index set ---" << std::endl;
            std::cout << "       Threshbound:         " << threshbound <<  std::endl;
            std::cout << "       ||P_{Lambda}r ||_2:  " << std::sqrt(P_Lambda_ResNE_square) << std::endl;
            std::cout << "       alpha*Residual_NE:   " << awgm_params.alpha*resNE_norm << std::endl;
        }

		// Add buckets to dummy vector: p
		for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
			P_Lambda_ResNE_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
			r_bucket.addBucketToCoefficients(p,i);

			if(awgm_params.verbose){
				std::cerr << "          Bucket " << i << ": L2-norm " << r_bucket.bucket_ell2norms[i] << std::endl;
			}

			if (P_Lambda_ResNE_square >= awgm_params.alpha*resNE_norm*awgm_params.alpha*resNE_norm) {
				break;
			}
		}

		// Replace indizes from last iteration
		for(auto& entry : u){
			resNE[entry.first] = 0.;
		}

		// New LambdaTrial: Add new indizes to u, complete to multitree
        IndexSet<Index2D> LambdaTrialDiff;
        for(auto& coeff : p){
			if(u.find(coeff.first) == u.end()){
				completeMultiTree(trialbasis, coeff.first, u, LambdaTrialDiff, 0, true);
			}
		}

        // Compute new u_leafs
        u_leafs.clear();
        FillWithZeros(LambdaTrialDiff,u_leafs);

        // New LambdaTest: stable expansion of LambdaTrial (on dummy vector Ap)
        Ap.clear();
        switch(awgm_params.stable_exp_u){
			case FullExpansion:
				getStableExpansion(trialbasis, testbasis, u, Ap);
				break;
			case WoMixedHWExpansion:
				getStableExpansion_woMixedHW(trialbasis, testbasis, u, Ap);
				break;
			case OnlyTemporalHWExpansion:
				getStableExpansion_onlyTemporalHW(trialbasis, testbasis, u, Ap);
				break;
			default:
				std::cerr << "Stable Expansion type doesn't exist!" << std::endl;
				break;
        }

        if(awgm_params.verbose){
            std::cout << "       LambdaTrial:  raw extension  " <<  std::setw(8) << std::right <<  p.size() << std::endl;
            std::cout << "                     multitree ext. " <<  std::setw(8) <<  LambdaTrialDiff.size() << std::endl;
            std::cout << "                     total size     " <<  std::setw(8) <<  u.size() << std::endl;
            std::cout << "       LambdaTest:   total size     " <<  std::setw(8) <<  Ap.size() << std::left << std::endl;
            if(awgm_params.verbose_extra){
            	std::cout << LambdaTrialDiff;
            }
            std::cout << std::endl;
        }
    }

    std::cerr << "AWGM reached maximal iteration number " << awgm_params.max_its << ": "
    		  << "Residual NE = " << resNE_norm << " "
    		  << ", Residual Au-f = " << res_norm << " "
              << ", awgm_tol = " << awgm_params.tol << std::endl << std::endl;

    if(awgm_params.print_info){
    	if(awgm_params.verbose){
            std::cout << "=====>  Writing AWGM CGLS information to file " << std::endl << std::endl;
    	}
    	awgm_info.print(awgm_params.info_filename.c_str());
    }
    if(awgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
    	if(awgm_params.verbose){
            std::cout << "=====>  Plotting AWGM CGLS solution to file " << std::endl << std::endl;
    	}
        plot2D<T,TrialBasis,TrialPrec>(trialbasis, u, trialPrec, exact_sol, 0., 1., 0., 1., 0.01, awgm_params.plot_filename.c_str());
    }
    
    if(awgm_params.clear_solver){
        Op.clear();
        OpTransp.clear();
        F.clear();
    }
    
    return;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
set_sol(sol_fct_2d _sol)
{
	exact_sol = _sol;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
AWGM_PG_Parameters&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
access_params()
{
	return awgm_params;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
reset_info()
{
	return awgm_info.reset();
}

} // namespace lawa
