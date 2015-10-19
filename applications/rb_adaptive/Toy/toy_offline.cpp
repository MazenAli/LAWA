#include "toy.h"


struct precon
{
    T
    operator() (const Index2D&)
    {
        return 1.;
    }
};

T
zero(T, T)
{
    return 0.;
}


int main () {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//

    int d   = 4;
    int j0  = 2;
    cout << "Wavelets: d = " << d <<  ", j0 = " << j0 << endl << endl; 


    /// Basis initialization
    _Basis     basis(d, j0);
    Basis2D    basis2d(basis, basis);


/*
    for (int i = 1; i<= 6; ++i) {
        DataType        u;

        std::string     snapshot;
        std::string     plot;
        snapshot = "training_data_toy/snap/snap_";
        plot     = "plot_";
        snapshot += std::to_string(i);
        plot     += std::to_string(i);
        snapshot += ".txt";

        std::cout << "Reading snapshot from file " << snapshot << std::endl;
        std::cout << "Printing plot to file " << plot << std::endl;
        readCoeffVector2D(u, snapshot.c_str());
        precon _p;
        plot2D(basis2d, u, _p,
                                 zero,
                                 0., 1.,
                                 0., 1.,
                                 1e-02,
                                 plot.c_str());
    }
    exit(0);
*/


    /// Initialization of operators
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);

    // Bilinear Forms
    Identity1D 		    IdentityBil(basis);
    Laplace1D 	        LaplaceBil(basis);

    RefIdentity1D 		    RefIdentityBil(basis.refinementbasis);
    RefLaplace1D 	        RefLaplaceBil(basis.refinementbasis);

    /// Initialization of local operator
    LOp_Id1D              lOp_Id1D      (basis, basis, RefIdentityBil, IdentityBil);
    LOp_Lapl1D            lOp_Lapl1D    (basis, basis, RefLaplaceBil,  LaplaceBil);

    LOp_Id_Id_2D			localIdentityIdentityOp2D		(lOp_Id1D, 		lOp_Id1D);
    LOp_Id_Lapl_2D			localIdentityLaplaceOp2D		(lOp_Id1D, 		lOp_Lapl1D);
    LOp_Lapl_Id_2D			localLaplaceIdentityOp2D		(lOp_Lapl1D, 		lOp_Id1D);

    localIdentityIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    localLaplaceIdentityOp2D.setJ(9);

    /// Initialization of preconditioner
    Prec2D prec(basis2d);
    NoPrec2D noPrec;

    /// Initialization of rhs

    /// Right Hand Side:
    ///     No Singular Supports in both dimensions
    DenseVectorT sing_support;
    ///      Forcing Functions
    Function2D<T> F1_Fct(f1, sing_support, sing_support);
    Function2D<T> F2_Fct(f2, sing_support, sing_support);
    Function2D<T> F3_Fct(f3, sing_support, sing_support);
    Function2D<T> F4_Fct(f4, sing_support, sing_support);
    Function2D<T> F5_Fct(f5, sing_support, sing_support);
    Function2D<T> F6_Fct(f6, sing_support, sing_support);
    Function2D<T> F7_Fct(f7, sing_support, sing_support);
    Function2D<T> F8_Fct(f8, sing_support, sing_support);
    Function2D<T> F9_Fct(f9, sing_support, sing_support);

    RhsIntegral2D			rhs_1(basis2d, F1_Fct,  100);
    RhsIntegral2D			rhs_2(basis2d, F2_Fct,  100);
    RhsIntegral2D			rhs_3(basis2d, F3_Fct,  100);
    RhsIntegral2D			rhs_4(basis2d, F4_Fct,  100);
    RhsIntegral2D			rhs_5(basis2d, F5_Fct,  100);
    RhsIntegral2D			rhs_6(basis2d, F6_Fct,  100);
    RhsIntegral2D			rhs_7(basis2d, F7_Fct,  100);
    RhsIntegral2D			rhs_8(basis2d, F8_Fct,  100);
    RhsIntegral2D			rhs_9(basis2d, F9_Fct,  100);
    Rhs           			F1(rhs_1,noPrec);
    Rhs           			F2(rhs_2,noPrec);
    Rhs           			F3(rhs_3,noPrec);
    Rhs           			F4(rhs_4,noPrec);
    Rhs           			F5(rhs_5,noPrec);
    Rhs           			F6(rhs_6,noPrec);
    Rhs           			F7(rhs_7,noPrec);
    Rhs           			F8(rhs_8,noPrec);
    Rhs           			F9(rhs_9,noPrec);

	//===============================================================//
	//===============  RB SETUP =====================================//
	//===============================================================//

    // Affine Decompositions:
    // 	Left Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct>	lhs_theta_fcts;
    lhs_theta_fcts.push_back(no_theta);
    ThetaStructure<ParamType> lhs_theta(lhs_theta_fcts);

    vector<AbstractLocalOperator2D<T>* > lhs_ops;
    lhs_ops.push_back(&localLaplaceIdentityOp2D);
    lhs_ops.push_back(&localIdentityLaplaceOp2D);
    Flex_COp_2D                       A(lhs_ops);
    vector<Flex_COp_2D*>                 ops_vec;
    ops_vec.push_back(&A);

    Affine_Op_2D affine_lhs(lhs_theta, ops_vec);

    H1_InnProd_2D innprod(localIdentityIdentityOp2D, localLaplaceIdentityOp2D, localIdentityLaplaceOp2D);

    // Right Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct> rhs_theta_fcts;
    rhs_theta_fcts.push_back(theta_chi_1);
    rhs_theta_fcts.push_back(theta_chi_2);
    rhs_theta_fcts.push_back(theta_chi_3);
    rhs_theta_fcts.push_back(theta_chi_4);
    rhs_theta_fcts.push_back(theta_chi_5);
    rhs_theta_fcts.push_back(theta_chi_6);
    rhs_theta_fcts.push_back(theta_chi_7);
    rhs_theta_fcts.push_back(theta_chi_8);
    rhs_theta_fcts.push_back(theta_chi_9);
    ThetaStructure<ParamType> rhs_theta(rhs_theta_fcts);
    vector<Rhs*> rhs_fcts;
    rhs_fcts.push_back(&F1);
    rhs_fcts.push_back(&F2);
    rhs_fcts.push_back(&F3);
    rhs_fcts.push_back(&F4);
    rhs_fcts.push_back(&F5);
    rhs_fcts.push_back(&F6);
    rhs_fcts.push_back(&F7);
    rhs_fcts.push_back(&F8);
    rhs_fcts.push_back(&F9);

    Affine_Rhs_2D affine_rhs(rhs_theta, rhs_fcts);
    RieszF_Rhs_2D rieszF_rhs(rhs_fcts);
    RieszA_Rhs_2D rieszA_rhs(ops_vec);

    // Right Hand Sides for direct Riesz Representor
    AffineA_Rhs_2D	 	affineA_rhs(lhs_theta, lhs_ops);
    RieszRes_Rhs_2D		rieszRes_rhs(affineA_rhs, affine_rhs);

	//===============================================================//
	//===============  AWGM =========================================//
	//===============================================================//

    IS_Parameters is_parameters;
    AWGM_Parameters awgm_truth_parameters, awgm_riesz_f_parameters, awgm_riesz_a_parameters, awgm_riesz_res_parameters;
    awgm_truth_parameters.max_its = 1e+03;

    is_parameters.adaptive_tol = true;
    is_parameters.init_tol = 0.0001;
    is_parameters.res_reduction = 0.01;

    //----------- Solver ---------------- //

    T gamma = 0.2;
    IndexSet<Index2D> Lambda;
    getSparseGridIndexSet(basis2d, Lambda, 1, 0, gamma);

    MT_AWGM_Truth awgm_u(basis2d, affine_lhs, affine_rhs, prec, awgm_truth_parameters, is_parameters);
    awgm_u.set_sol(dummy);
    awgm_u.awgm_params.tol = 1e-04;
    awgm_u.awgm_params.max_basissize = 1e+06;
    awgm_u.set_initial_indexset(Lambda);

    MT_AWGM_Riesz_F awgm_rieszF(basis2d, innprod, rieszF_rhs, prec, awgm_riesz_f_parameters, is_parameters);
    awgm_rieszF.set_sol(dummy);
    awgm_rieszF.set_initial_indexset(Lambda);
    awgm_rieszF.awgm_params.tol = 5e-04;
    awgm_rieszF.awgm_params.info_filename = "stempel_rieszF_conv_info.txt";

    MT_AWGM_Riesz_A awgm_rieszA(basis2d, innprod, rieszA_rhs, prec, awgm_riesz_a_parameters, is_parameters);
    awgm_rieszA.set_sol(dummy);
    awgm_rieszA.set_initial_indexset(Lambda);
    awgm_rieszA.awgm_params.tol = 5e-04;
    awgm_rieszA.awgm_params.info_filename = "stempel_rieszA_conv_info.txt";

    MT_AWGM_Riesz_Res awgm_rieszRes(basis2d, innprod, rieszRes_rhs, prec, awgm_riesz_res_parameters, is_parameters);
    awgm_rieszRes.set_sol(dummy);
    awgm_rieszRes.set_initial_indexset(Lambda);
    awgm_rieszRes.awgm_params.tol = 5e-04;
    awgm_rieszRes.awgm_params.print_info = false;

    MTTruthSolver rb_truth(awgm_u, awgm_rieszF, awgm_rieszA, &awgm_rieszRes, innprod, affine_lhs, rieszF_rhs);

    //----------- RB System ---------------- //

    RB_Model rb_system(lhs_theta, rhs_theta);

    rb_system.rb_params.ref_param = {{1.}};
    rb_system.rb_params.call = call_cg;

    //----------- RB Base ---------------- //

    RB_BaseModel rb_base(rb_system, rb_truth, Lambda);

    ParamType mu_min = {{0}};
    ParamType mu_max = {{1}};

    rb_base.greedy_params.training_type = weak;
    rb_base.greedy_params.tol = 1e-05;
    rb_base.greedy_params.min_param = mu_min;
    rb_base.greedy_params.max_param = mu_max;
    rb_base.greedy_params.Nmax =    6;
    rb_base.greedy_params.nb_training_params = {{nb_stempel}};
    rb_base.greedy_params.log_scaling = {{0}};
    rb_base.greedy_params.print_file = "toy_greedy_info.txt";
    rb_base.greedy_params.trainingdata_folder = "training_data_toy";
    rb_base.greedy_params.print_paramset = true;
    rb_base.greedy_params.erase_snapshot_params = false;
    rb_base.greedy_params.orthonormalize_bfs = true;
    rb_base.greedy_params.tighten_tol = true;
    rb_base.greedy_params.snapshot_tol_red_crit = repeated_param;
    rb_base.greedy_params.tighten_tol_rieszA = false;
    rb_base.greedy_params.tighten_tol_rieszF = false;
    rb_base.greedy_params.tighten_tol_reduction = 0.1,
    rb_base.greedy_params.update_snapshot = true;
    rb_base.greedy_params.update_rieszF = false;
    rb_base.greedy_params.update_rieszA = false;
    rb_base.greedy_params.test_estimator_equivalence = false;
    rb_base.greedy_params.write_direct_representors = false;
    rb_base.greedy_params.riesz_constant_X = 5.6;
    rb_base.greedy_params.riesz_constant_Y = 2;

    cout << "Parameters Training: " << std::endl << std::endl;
    rb_base.greedy_params.print();
    rb_system.rb_params.print();

    cout << "Parameters Truth Solver: " << std::endl << std::endl;
    awgm_u.awgm_params.print();

    awgm_u.is_params.print();

    cout << "Parameters Riesz Solver F : " << std::endl << std::endl;
    awgm_rieszF.awgm_params.print();

    cout << "Parameters Riesz Solver A : " << std::endl << std::endl;
    awgm_rieszA.awgm_params.print();

    cout << "Parameters Riesz Solver Res : " << std::endl << std::endl;
    awgm_rieszRes.awgm_params.print();

    rb_base.init();
    rb_base.train_Greedy();

    return 0;
}

