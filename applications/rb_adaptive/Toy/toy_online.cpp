#include "toy.h"
#include <fstream>
#include <sstream>
#include <ctime>


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


void
read_paramfile(string paramfilename, vector<ParamType>& v){
    T mu1;
    ifstream paramfile(paramfilename);
    if(paramfile.is_open()){
    	while(paramfile >> mu1){
        	ParamType mu = {{mu1}};
        	v.push_back(mu);
    	}
    	paramfile.close();
    }
    else{
    	cerr << "Couldn't open " << paramfilename << " for reading!" << endl;
    }
}

int main (int argc, char* argv[]) {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//
	
    int d   = 4;
    int j0  = 2;

    if(argc < 5 || argc > 7){
    	cerr << "Usage: " << argv[0] << " offline_data_folder output_name readSols(0/1) updateSols(0/1) (InputSolutionFolder) (OutputSolutionFolder)" << endl;
    	exit(1);
     }

    string offline_folder = argv[1];
    string output = argv[2];
    int readSols = atoi(argv[3]);
    int updateSols = atoi(argv[4]);
    string solfolder_in, solfolder_out;
    if(argc >= 6){
        solfolder_in = argv[5];
    }
    if(argc == 7){
        solfolder_out = argv[6];
    }

    /// Basis initialization
    _Basis     basis(d,j0);
    Basis2D basis2d(basis,basis);
    
    /// Initialization of operators
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);

    // Bilinear Forms
    Identity1D  		    IdentityBil(basis);
    Laplace1D   	        LaplaceBil(basis);
    RefIdentity1D    		RefIdentityBil(basis.refinementbasis);
    RefLaplace1D 	        RefLaplaceBil(basis.refinementbasis);

    /// Initialization of local operator
    LOp_Id1D              lOp_Id1D      (basis, basis, RefIdentityBil,   IdentityBil);
    LOp_Lapl1D            lOp_Lapl1D    (basis, basis, RefLaplaceBil,    LaplaceBil);

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
    Function2D<T> F10_Fct(f10, sing_support, sing_support);
    Function2D<T> F11_Fct(f11, sing_support, sing_support);
    Function2D<T> F12_Fct(f12, sing_support, sing_support);
    Function2D<T> F13_Fct(f13, sing_support, sing_support);
    Function2D<T> F14_Fct(f14, sing_support, sing_support);
    Function2D<T> F15_Fct(f15, sing_support, sing_support);
    Function2D<T> F16_Fct(f16, sing_support, sing_support);
    Function2D<T> F17_Fct(f17, sing_support, sing_support);
    Function2D<T> F18_Fct(f18, sing_support, sing_support);
    Function2D<T> F19_Fct(f19, sing_support, sing_support);
    Function2D<T> F20_Fct(f20, sing_support, sing_support);
    Function2D<T> F21_Fct(f21, sing_support, sing_support);
    Function2D<T> F22_Fct(f22, sing_support, sing_support);
    Function2D<T> F23_Fct(f23, sing_support, sing_support);
    Function2D<T> F24_Fct(f24, sing_support, sing_support);
    Function2D<T> F25_Fct(f25, sing_support, sing_support);
    Function2D<T> F26_Fct(f26, sing_support, sing_support);
    Function2D<T> F27_Fct(f27, sing_support, sing_support);
    Function2D<T> F28_Fct(f28, sing_support, sing_support);
    Function2D<T> F29_Fct(f29, sing_support, sing_support);
    Function2D<T> F30_Fct(f30, sing_support, sing_support);
    Function2D<T> F31_Fct(f31, sing_support, sing_support);
    Function2D<T> F32_Fct(f32, sing_support, sing_support);
    Function2D<T> F33_Fct(f33, sing_support, sing_support);
    Function2D<T> F34_Fct(f34, sing_support, sing_support);
    Function2D<T> F35_Fct(f35, sing_support, sing_support);
    Function2D<T> F36_Fct(f36, sing_support, sing_support);
    Function2D<T> F37_Fct(f37, sing_support, sing_support);
    Function2D<T> F38_Fct(f38, sing_support, sing_support);
    Function2D<T> F39_Fct(f39, sing_support, sing_support);
    Function2D<T> F40_Fct(f40, sing_support, sing_support);
    Function2D<T> F41_Fct(f41, sing_support, sing_support);
    Function2D<T> F42_Fct(f42, sing_support, sing_support);
    Function2D<T> F43_Fct(f43, sing_support, sing_support);
    Function2D<T> F44_Fct(f44, sing_support, sing_support);
    Function2D<T> F45_Fct(f45, sing_support, sing_support);
    Function2D<T> F46_Fct(f46, sing_support, sing_support);
    Function2D<T> F47_Fct(f47, sing_support, sing_support);
    Function2D<T> F48_Fct(f48, sing_support, sing_support);
    Function2D<T> F49_Fct(f49, sing_support, sing_support);


    ///     Peaks: points and corresponding coefficients
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    RhsIntegral2D			rhs_1(basis2d, F1_Fct,  100);
    RhsIntegral2D			rhs_2(basis2d, F2_Fct,  100);
    RhsIntegral2D			rhs_3(basis2d, F3_Fct,  100);
    RhsIntegral2D			rhs_4(basis2d, F4_Fct,  100);
    RhsIntegral2D			rhs_5(basis2d, F5_Fct,  100);
    RhsIntegral2D			rhs_6(basis2d, F6_Fct,  100);
    RhsIntegral2D			rhs_7(basis2d, F7_Fct,  100);
    RhsIntegral2D			rhs_8(basis2d, F8_Fct,  100);
    RhsIntegral2D			rhs_9(basis2d, F9_Fct,  100);
    RhsIntegral2D			rhs_10(basis2d, F10_Fct,  100);
    RhsIntegral2D			rhs_11(basis2d, F11_Fct,  100);
    RhsIntegral2D			rhs_12(basis2d, F12_Fct,  100);
    RhsIntegral2D			rhs_13(basis2d, F13_Fct,  100);
    RhsIntegral2D			rhs_14(basis2d, F14_Fct,  100);
    RhsIntegral2D			rhs_15(basis2d, F15_Fct,  100);
    RhsIntegral2D			rhs_16(basis2d, F16_Fct,  100);
    RhsIntegral2D			rhs_17(basis2d, F17_Fct,  100);
    RhsIntegral2D			rhs_18(basis2d, F18_Fct,  100);
    RhsIntegral2D			rhs_19(basis2d, F19_Fct,  100);
    RhsIntegral2D			rhs_20(basis2d, F20_Fct,  100);
    RhsIntegral2D			rhs_21(basis2d, F21_Fct,  100);
    RhsIntegral2D			rhs_22(basis2d, F22_Fct,  100);
    RhsIntegral2D			rhs_23(basis2d, F23_Fct,  100);
    RhsIntegral2D			rhs_24(basis2d, F24_Fct,  100);
    RhsIntegral2D			rhs_25(basis2d, F25_Fct,  100);
    RhsIntegral2D			rhs_26(basis2d, F26_Fct,  100);
    RhsIntegral2D			rhs_27(basis2d, F27_Fct,  100);
    RhsIntegral2D			rhs_28(basis2d, F28_Fct,  100);
    RhsIntegral2D			rhs_29(basis2d, F29_Fct,  100);
    RhsIntegral2D			rhs_30(basis2d, F30_Fct,  100);
    RhsIntegral2D			rhs_31(basis2d, F31_Fct,  100);
    RhsIntegral2D			rhs_32(basis2d, F32_Fct,  100);
    RhsIntegral2D			rhs_33(basis2d, F33_Fct,  100);
    RhsIntegral2D			rhs_34(basis2d, F34_Fct,  100);
    RhsIntegral2D			rhs_35(basis2d, F35_Fct,  100);
    RhsIntegral2D			rhs_36(basis2d, F36_Fct,  100);
    RhsIntegral2D			rhs_37(basis2d, F37_Fct,  100);
    RhsIntegral2D			rhs_38(basis2d, F38_Fct,  100);
    RhsIntegral2D			rhs_39(basis2d, F39_Fct,  100);
    RhsIntegral2D			rhs_40(basis2d, F40_Fct,  100);
    RhsIntegral2D			rhs_41(basis2d, F41_Fct,  100);
    RhsIntegral2D			rhs_42(basis2d, F42_Fct,  100);
    RhsIntegral2D			rhs_43(basis2d, F43_Fct,  100);
    RhsIntegral2D			rhs_44(basis2d, F44_Fct,  100);
    RhsIntegral2D			rhs_45(basis2d, F45_Fct,  100);
    RhsIntegral2D			rhs_46(basis2d, F46_Fct,  100);
    RhsIntegral2D			rhs_47(basis2d, F47_Fct,  100);
    RhsIntegral2D			rhs_48(basis2d, F48_Fct,  100);
    RhsIntegral2D			rhs_49(basis2d, F49_Fct,  100);

    Rhs           			F1(rhs_1,noPrec);
    Rhs           			F2(rhs_2,noPrec);
    Rhs           			F3(rhs_3,noPrec);
    Rhs           			F4(rhs_4,noPrec);
    Rhs           			F5(rhs_5,noPrec);
    Rhs           			F6(rhs_6,noPrec);
    Rhs           			F7(rhs_7,noPrec);
    Rhs           			F8(rhs_8,noPrec);
    Rhs           			F9(rhs_9,noPrec);
    Rhs           			F10(rhs_10,noPrec);
    Rhs           			F11(rhs_11,noPrec);
    Rhs           			F12(rhs_12,noPrec);
    Rhs           			F13(rhs_13,noPrec);
    Rhs           			F14(rhs_14,noPrec);
    Rhs           			F15(rhs_15,noPrec);
    Rhs           			F16(rhs_16,noPrec);
    Rhs           			F17(rhs_17,noPrec);
    Rhs           			F18(rhs_18,noPrec);
    Rhs           			F19(rhs_19,noPrec);
    Rhs           			F20(rhs_20,noPrec);
    Rhs           			F21(rhs_21,noPrec);
    Rhs           			F22(rhs_22,noPrec);
    Rhs           			F23(rhs_23,noPrec);
    Rhs           			F24(rhs_24,noPrec);
    Rhs           			F25(rhs_25,noPrec);
    Rhs           			F26(rhs_26,noPrec);
    Rhs           			F27(rhs_27,noPrec);
    Rhs           			F28(rhs_28,noPrec);
    Rhs           			F29(rhs_29,noPrec);
    Rhs           			F30(rhs_30,noPrec);
    Rhs           			F31(rhs_31,noPrec);
    Rhs           			F32(rhs_32,noPrec);
    Rhs           			F33(rhs_33,noPrec);
    Rhs           			F34(rhs_34,noPrec);
    Rhs           			F35(rhs_35,noPrec);
    Rhs           			F36(rhs_36,noPrec);
    Rhs           			F37(rhs_37,noPrec);
    Rhs           			F38(rhs_38,noPrec);
    Rhs           			F39(rhs_39,noPrec);
    Rhs           			F40(rhs_40,noPrec);
    Rhs           			F41(rhs_41,noPrec);
    Rhs           			F42(rhs_42,noPrec);
    Rhs           			F43(rhs_43,noPrec);
    Rhs           			F44(rhs_44,noPrec);
    Rhs           			F45(rhs_45,noPrec);
    Rhs           			F46(rhs_46,noPrec);
    Rhs           			F47(rhs_47,noPrec);
    Rhs           			F48(rhs_48,noPrec);
    Rhs           			F49(rhs_49,noPrec);



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
    rhs_theta_fcts.push_back(theta_chi_10);
    rhs_theta_fcts.push_back(theta_chi_11);
    rhs_theta_fcts.push_back(theta_chi_12);
    rhs_theta_fcts.push_back(theta_chi_13);
    rhs_theta_fcts.push_back(theta_chi_14);
    rhs_theta_fcts.push_back(theta_chi_15);
    rhs_theta_fcts.push_back(theta_chi_16);
    rhs_theta_fcts.push_back(theta_chi_17);
    rhs_theta_fcts.push_back(theta_chi_18);
    rhs_theta_fcts.push_back(theta_chi_19);
    rhs_theta_fcts.push_back(theta_chi_20);
    rhs_theta_fcts.push_back(theta_chi_21);
    rhs_theta_fcts.push_back(theta_chi_22);
    rhs_theta_fcts.push_back(theta_chi_23);
    rhs_theta_fcts.push_back(theta_chi_24);
    rhs_theta_fcts.push_back(theta_chi_25);
    rhs_theta_fcts.push_back(theta_chi_26);
    rhs_theta_fcts.push_back(theta_chi_27);
    rhs_theta_fcts.push_back(theta_chi_28);
    rhs_theta_fcts.push_back(theta_chi_29);
    rhs_theta_fcts.push_back(theta_chi_30);
    rhs_theta_fcts.push_back(theta_chi_31);
    rhs_theta_fcts.push_back(theta_chi_32);
    rhs_theta_fcts.push_back(theta_chi_33);
    rhs_theta_fcts.push_back(theta_chi_34);
    rhs_theta_fcts.push_back(theta_chi_35);
    rhs_theta_fcts.push_back(theta_chi_36);
    rhs_theta_fcts.push_back(theta_chi_37);
    rhs_theta_fcts.push_back(theta_chi_38);
    rhs_theta_fcts.push_back(theta_chi_39);
    rhs_theta_fcts.push_back(theta_chi_40);
    rhs_theta_fcts.push_back(theta_chi_41);
    rhs_theta_fcts.push_back(theta_chi_42);
    rhs_theta_fcts.push_back(theta_chi_43);
    rhs_theta_fcts.push_back(theta_chi_44);
    rhs_theta_fcts.push_back(theta_chi_45);
    rhs_theta_fcts.push_back(theta_chi_46);
    rhs_theta_fcts.push_back(theta_chi_47);
    rhs_theta_fcts.push_back(theta_chi_48);
    rhs_theta_fcts.push_back(theta_chi_49);

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
    rhs_fcts.push_back(&F10);
    rhs_fcts.push_back(&F11);
    rhs_fcts.push_back(&F12);
    rhs_fcts.push_back(&F13);
    rhs_fcts.push_back(&F14);
    rhs_fcts.push_back(&F15);
    rhs_fcts.push_back(&F16);
    rhs_fcts.push_back(&F17);
    rhs_fcts.push_back(&F18);
    rhs_fcts.push_back(&F19);
    rhs_fcts.push_back(&F20);
    rhs_fcts.push_back(&F21);
    rhs_fcts.push_back(&F22);
    rhs_fcts.push_back(&F23);
    rhs_fcts.push_back(&F24);
    rhs_fcts.push_back(&F25);
    rhs_fcts.push_back(&F26);
    rhs_fcts.push_back(&F27);
    rhs_fcts.push_back(&F28);
    rhs_fcts.push_back(&F29);
    rhs_fcts.push_back(&F30);
    rhs_fcts.push_back(&F31);
    rhs_fcts.push_back(&F32);
    rhs_fcts.push_back(&F33);
    rhs_fcts.push_back(&F34);
    rhs_fcts.push_back(&F35);
    rhs_fcts.push_back(&F36);
    rhs_fcts.push_back(&F37);
    rhs_fcts.push_back(&F38);
    rhs_fcts.push_back(&F39);
    rhs_fcts.push_back(&F40);
    rhs_fcts.push_back(&F41);
    rhs_fcts.push_back(&F42);
    rhs_fcts.push_back(&F43);
    rhs_fcts.push_back(&F44);
    rhs_fcts.push_back(&F45);
    rhs_fcts.push_back(&F46);
    rhs_fcts.push_back(&F47);
    rhs_fcts.push_back(&F48);
    rhs_fcts.push_back(&F49);

    Affine_Rhs_2D affine_rhs(rhs_theta, rhs_fcts);
    RieszF_Rhs_2D rieszF_rhs(rhs_fcts);
    RieszA_Rhs_2D rieszA_rhs(ops_vec);


	//===============================================================//
	//===============  AWGM =========================================//
	//===============================================================//


    IS_Parameters is_parameters;
    AWGM_Parameters awgm_truth_parameters, awgm_riesz_f_parameters, awgm_riesz_a_parameters;
    awgm_truth_parameters.max_its = 1000;
    
    is_parameters.adaptive_tol = true;
    is_parameters.init_tol = 0.0001;
    is_parameters.res_reduction = 0.01;

    //----------- Solver ---------------- //

    T gamma = 0.2;
    IndexSet<Index2D> Lambda;
    getSparseGridIndexSet(basis2d, Lambda, 1,0,gamma);

    MT_AWGM_Truth awgm_u(basis2d, affine_lhs, affine_rhs, prec, awgm_truth_parameters, is_parameters);
    awgm_u.set_sol(dummy);
    awgm_u.awgm_params.tol = 3e-06;
    awgm_u.awgm_params.max_basissize = 3e+06;
    awgm_u.set_initial_indexset(Lambda);

    MT_AWGM_Riesz_F awgm_rieszF(basis2d, innprod, rieszF_rhs, prec, awgm_riesz_f_parameters, is_parameters);
    awgm_rieszF.set_sol(dummy);
    awgm_rieszF.set_initial_indexset(Lambda);
    awgm_rieszF.awgm_params.tol = 5e-03;
    awgm_rieszF.awgm_params.info_filename = "stempel_rieszF_conv_info.txt";

    MT_AWGM_Riesz_A awgm_rieszA(basis2d, innprod, rieszA_rhs, prec, awgm_riesz_a_parameters, is_parameters);
    awgm_rieszA.set_sol(dummy);
    awgm_rieszA.set_initial_indexset(Lambda);
    awgm_rieszA.awgm_params.tol = 5e-03;
    awgm_rieszA.awgm_params.info_filename = "stempel_rieszA_conv_info.txt";

    MTTruthSolver rb_truth(awgm_u, awgm_rieszF, awgm_rieszA, innprod, affine_lhs, rieszF_rhs);

    //----------- RB System ---------------- //

    RB_Model rb_system(lhs_theta, rhs_theta);

    rb_system.rb_params.ref_param = {{1.}};
    rb_system.rb_params.call = call_cg;

    //----------- RB Base ---------------- //

    RB_BaseModel rb_base(rb_system, rb_truth, Lambda);


    ParamType mu_min = {{0}};
    ParamType mu_max = {{2}};

    rb_base.greedy_params.tol = 1e-05;
    rb_base.greedy_params.min_param = mu_min;
    rb_base.greedy_params.max_param = mu_max;
    rb_base.greedy_params.Nmax = 	6;
    rb_base.greedy_params.nb_training_params = {{nb_stempel}};
    rb_base.greedy_params.log_scaling = {{0}};
    rb_base.greedy_params.print_file = "toy_greedy_info.txt";
    rb_base.greedy_params.trainingdata_folder = "training_data_toy";
    rb_base.greedy_params.print_paramset = true;
    rb_base.greedy_params.erase_snapshot_params = false;
    rb_base.greedy_params.orthonormalize_bfs = false;
    rb_base.greedy_params.tighten_tol = false;
    rb_base.greedy_params.tighten_tol_rieszA = true;
    rb_base.greedy_params.tighten_tol_rieszF = true;
    rb_base.greedy_params.update_snapshot = true;
    rb_base.greedy_params.update_rieszF = true;
    rb_base.greedy_params.update_rieszA = false;
    rb_base.greedy_params.test_estimator_equivalence = true;
	
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

    // Read RB Data
    rb_system.read_rb_data("training_data_toy");
    rb_base.read_basisfunctions("training_data_toy/bf");

    // Read Test Parameters
    std::vector<ParamType> Xi_test;
    read_paramfile("Xitest.txt", Xi_test);


    cout << "===============================" << endl;
    cout << "Test Parameters: " << endl;
    for(auto& mu : Xi_test){
    	cout << mu[0] << endl;
    }
    cout << "===============================" << endl;


    // Test Reduced Solutions
    std::map<ParamType, std::vector<T> > errs, errbounds;
    int Nmax = rb_system.RB_inner_product.numRows();

    DenseVectorT av_err(Nmax);
    DenseVectorT av_errest(Nmax);


    // Reference solutions
    std::vector<T (*)(T)>   refvec;
    std::vector<T (*)(T)>   refvecdx;
    refvec.push_back(u1);
    refvec.push_back(u2);
    refvec.push_back(u3);
    refvec.push_back(u4);
    refvec.push_back(u5);
    refvec.push_back(u6);
    refvec.push_back(u7);
    refvecdx.push_back(ux1);
    refvecdx.push_back(ux2);
    refvecdx.push_back(ux3);
    refvecdx.push_back(ux4);
    refvecdx.push_back(ux5);
    refvecdx.push_back(ux6);
    refvecdx.push_back(ux7);


    int     count     = 0;
    int     num_old   = 3;
    int     num_new   = 4;
    int     num_tot   = num_old+num_new;
    double  duration  = 0.;
	// For all parameters:
    for(auto& mu : Xi_test){
        ++count;
        cout << "Mu = [" << mu[0]  << "]" << endl;

    	stringstream ufilename_in, ufilename_out;
    	if(argc >= 6){
    		ufilename_in << solfolder_in << "/";
    	}
    	if(argc == 7){
            ufilename_out << solfolder_out << "/";
    	}
    	ufilename_in << "u_" << mu[0]  << ".txt";
    	ufilename_out << "u_" << mu[0] << ".txt";

    	std::vector<T> errs_mu, errbounds_mu;

    	T err, errbound;

        int numx, numy;
        if (count <= num_old*num_old) {
            numy = (count-1)/num_old;
            numx = (count-1)%num_old;
        } else if ((count-num_old*num_old)<= num_old*num_new) {
            int temp = count - num_old*num_old;
            numx = (temp-1)/num_new;
            numy = (temp-1)%num_new;
            numy += num_old;
        } else {
            int temp = count - (num_old*num_old+num_old*num_new);
            numx = (temp-1)/num_tot;
            numy = (temp-1)%num_tot;
            numx += num_old;
        }

        cout << "Corresponding point (x,y) is ("
             << numx+1 << "," << numy+1 << ")\n";

		for(size_t n = Nmax; n <= Nmax; ++n){

    		// Compute RB solution
    		DenseVectorT u_N = rb_system.get_rb_solution(n, mu);
            cout << "u_" << n << " = " << u_N << endl;
    		// Compute real error
    		DataType u_approx = rb_base.reconstruct_u_N(u_N, n);
            // Sequential ordering 1...49 of parameters assumed
    		SeparableFunction2D<T>  uref (refvec[numx],   sing_support,
                                          refvec[numy],   sing_support);
    		SeparableFunction2D<T>  uxref(refvecdx[numx], sing_support,
                                          refvec[numy],   sing_support);
    		SeparableFunction2D<T>  uyref(refvec[numx],   sing_support,
                                          refvecdx[numy], sing_support);
            //

            err = h1Error2D(basis2d, uref, uxref, uyref, u_approx);
            cout << "Err = " << err << endl;
    		errs_mu.push_back(err);
    		
    		av_err(n) += err;

            // Compute error estimator
            std::clock_t t;
            t = std::clock();
    		errbound = rb_system.get_errorbound(u_N, mu);
            t = std::clock() - t;
            std::cout << "\n***It took " << t << " clocks to compute***\n";
            duration += t / (double) CLOCKS_PER_SEC;
    		errbounds_mu.push_back(errbound);
    		av_errest(n) += errbound;
    		
		}

    	errs.insert(make_pair(mu, errs_mu));
    	errbounds.insert(make_pair(mu, errbounds_mu));
    }

    std::cout << "Total   computation time is " << duration << " secs\n";
    duration /= Nmax;
    std::cout << "Average computation time is " << duration << " secs\n";

    av_err *= 1./Xi_test.size();
	av_errest *= 1./Xi_test.size();

    ofstream errfile(output + "_OnlineTest_Errors.txt");
    if(errfile.is_open()){
    	errfile << "Mu Err(N=1) Err(N=2) ...." << endl;
    	for(auto& entry : errs){
    		errfile << entry.first[0];
    		for(auto& e : entry.second){
    			errfile << " " << e;
    		}
    		errfile << endl;
    	}
    	errfile.close();
    }
    else{
    	cerr << "Couldn't open Test_Errors.txt for writing!" << endl;
    }

    ofstream errestfile(output + "_OnlineTest_ErrorEstimators.txt");
    if(errestfile.is_open()){
    	errestfile << "Mu1 ErrEst(N=1) ErrEst(N=2) ...." << endl;
    	for(auto& entry : errbounds){
    		errestfile << entry.first[0];
    		for(auto& e : entry.second){
    			errestfile << " " << e;
    		}
    		errestfile << endl;
    	}
    	errestfile.close();
    }
    else{
    	cerr << "Couldn't open Test_Errors.txt for writing!" << endl;
    }

    
    ofstream av_file(output + "_OnlineTest_Averages.txt");
    if(av_file.is_open()){
        av_file << "# N avErr avErrEst" << endl;
        for(int i = 1; i <= Nmax; ++i){
            av_file << i << " " << i << " " << av_err(i) << " " << av_errest(i) << endl;
        }
        
        av_file.close();        
    }
    else{
        cerr << "Couldn't open file " << output << "_OnlineTest_Averages.txt" << " for writing! " << endl;
    }

    std::cout << "Everything exited fine\n";

    return 0;
}

