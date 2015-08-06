#include <iostream>
#include <iomanip>
#include <utility>
#include <array> 
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

// quick and dirty finite wavelet expansion implementation
template <typename T, typename Index, typename RHSINTEGRAL,
          typename Preconditioner>
class SimpleRHS : public RHS<T, Index, RHSINTEGRAL, Preconditioner>
{
    private:
        const IndexSet<Index>&      LambdaTest;
    public:
        SimpleRHS(const RHSINTEGRAL& rhsintegral, Preconditioner& P,
                  IndexSet<Index>& _LambdaTest):
                  RHS<T, Index, RHSINTEGRAL, Preconditioner>(rhsintegral, P),
                  LambdaTest(_LambdaTest){}

        T
        operator()(const Index& lambda) override
        {
            typename IndexSet<Index>::const_iterator
                                      found = LambdaTest.find(lambda);
            if (found == LambdaTest.end()) {
                return 0.;
            } else {
                return RHS<T, Index, RHSINTEGRAL, Preconditioner>
                           ::operator()(lambda);
            }
        }

        Coefficients<Lexicographical,T,Index>
        operator()(const IndexSet<Index>& Lambda) override
        {
            typedef typename IndexSet<Index>::const_iterator const_set_it;
            Coefficients<Lexicographical,T,Index> ret;
            for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
                T tmp = SimpleRHS<T,Index,RHSINTEGRAL,Preconditioner>::operator()(*lambda);
                ret[*lambda] = tmp;
            }
            return ret;
        }
};

template <typename LocalOperator1, typename LocalOperator2>
struct SimpleOp2D : LocalOperator2D<LocalOperator1, LocalOperator2>
{
    typedef typename LocalOperator1::T      T;
    const IndexSet<Index2D>&                Lambda; // must be a multi-tree!

    SimpleOp2D(LocalOperator1& _localop1, LocalOperator2& _localop2,
               IndexSet<Index2D>& _Lambda):
               LocalOperator2D<LocalOperator1, LocalOperator2>
               (_localop1, _localop2), Lambda(_Lambda){}

    void
    eval(const Coefficients<Lexicographical, T, Index2D>& v,
               Coefficients<Lexicographical, T, Index2D>& AAv) override
    {
        Coefficients<Lexicographical, T, Index2D>   temp = v;

        P(Lambda, temp);
        P(Lambda, AAv);
        LocalOperator2D<LocalOperator1, LocalOperator2>::eval(temp, AAv);
    }
};

// quick and dirty implementation of new error estimator

template <typename RB_Model, typename TruthModel, typename DataType,
          typename ParamType, typename _IndexSet>
class RB_Base_wavest : public RB_Base<RB_Model, TruthModel, DataType, ParamType>
{
    private:
        // sparse grid and thresh parameter
        const _IndexSet&    Lambda;
        double              ETA     = 1e-16;
        bool                recalcF = false;
        long                total   = 0;

    public:
        RB_Base_wavest(RB_Model& _rb_system, TruthModel& _rb_truth,
                       const _IndexSet& _Lambda)
                      :RB_Base<RB_Model, TruthModel, DataType, ParamType>
                       (_rb_system, _rb_truth),
                       Lambda(_Lambda){}


        void
        init()
        {
            this->greedy_params.tighten_tol_rieszA              = false;
            this->greedy_params.tighten_tol_rieszF              = false;
            this->greedy_params.update_rieszF                   = false;
            this->greedy_params.update_rieszA                   = false;
            this->greedy_params.test_estimator_equivalence      = false;
            this->greedy_params.tighten_estimator_accuracy      = false;
            std::cout << "Exact Lambda size is " << Lambda.size() << std::endl;
        }

    private:
        void
        calculate_Riesz_RHS_information(bool) override
        {
            std::cout << "Correct RHS being called\n";
            if (recalcF) {
                std::cerr << "Warning: recalculating RHS! This should not happen!\n";
            }
            calc_wav_RHS();
            recalcF = true;

            std::size_t Qf = this->rb_system.Q_f();
            this->rb_system.F_F_representor_norms.engine().resize((int) Qf, (int) Qf);
            for(std::size_t qf1 = 1; qf1 <= Qf; ++qf1) {
                for (std::size_t qf2 = qf1; qf2 <= Qf; ++qf2) {

                    this->rb_system.F_F_representor_norms(qf1,qf2) =
                          dotProduct(this->F_representors[qf1-1],
                                     this->F_representors[qf2-1]);

                    if(qf1 != qf2) {
                        this->rb_system.F_F_representor_norms(qf2,qf1) = this->rb_system.F_F_representor_norms(qf1,qf2);
                    }
                }
            }
        }


        void
        calculate_Riesz_LHS_information(const DataType& bf, bool) override
        {
            std::cout << "Correct LHS called\n";
            int N = this->rb_basisfunctions.size();
            std::cout << "Current N=" << N << std::endl;
            calc_wav_LHS(bf);

            // Update the Riesz Representor Norms A x A
            std::size_t Qa = this->rb_system.Q_a();
            for(int n1 = 0; n1 < N; ++n1) {
                typename RB_Model::FullColMatrixT A_n1_N(Qa, Qa);
                for(std::size_t qa1 = 1; qa1 <= Qa; ++qa1) {
                    for(std::size_t qa2 = qa1; qa2 <= Qa; ++qa2) {
                        A_n1_N(qa1, qa2) =
                        dotProduct(this->A_representors[n1][qa1-1],
                                   this->A_representors[N-1][qa2-1]);

                        if(qa1 != qa2){
                            if(n1 == N-1){
                                A_n1_N(qa2, qa1) = A_n1_N(qa1, qa2);
                            }
                            else{
                                A_n1_N(qa2, qa1) =
                                dotProduct(this->A_representors[n1][qa2-1],
                                           this->A_representors[N-1][qa1-1]);
                            }
                        }
                    }
                }
                if(n1 == N-1){
                    std::vector<typename RB_Model::FullColMatrixT> newvec;
                    newvec.push_back(A_n1_N);
                    this->rb_system.A_A_representor_norms.push_back(newvec);
                }
                else{
                    this->rb_system.A_A_representor_norms[n1].push_back(A_n1_N);
                }
            }

            // Update the Riesz Representor Norms A x F
            std::size_t Qf = this->rb_system.Q_f();
            typename RB_Model::FullColMatrixT A_F(Qa, Qf);
            for(std::size_t qa = 1; qa <= Qa; ++qa) {
                for(std::size_t qf = 1; qf <= Qf; ++qf) {
                    A_F(qa, qf) =
                    dotProduct(this->A_representors[N-1][qa-1],
                               this->F_representors[qf-1]);
                }
            }
            this->rb_system.A_F_representor_norms.push_back(A_F);
        }


        void
        calc_wav_RHS()
        {
            assert(this->F_representors.size()==0);

            auto& F = this->rb_truth.access_RieszSolver_F().get_rhs();
            auto& P = this->rb_truth.access_RieszSolver_F().get_testprec();

            std::size_t Qf = this->rb_system.Q_f();
            assert(Qf>0);
            this->F_representors.resize(Qf);

            for (std::size_t i=0; i<Qf; ++i) {
                F.set_active_comp(i);
                this->F_representors[i] = F(Lambda, P); // precondtioned
                THRESH_NoCopy(this->F_representors[i], ETA);
                this->total += this->F_representors[i].size();
            }
            std::cout << "Current total size of F is " << total << std::endl;
        }


        void
        calc_wav_LHS(const DataType& bf)
        {
            std::size_t             Qa = this->rb_system.Q_a();
            assert(Qa>0);
            std::vector<DataType>   temp(Qa);
            auto& A = this->rb_truth.access_RieszSolver_A().get_rhs().get_bilformvec();
            auto& P = this->rb_truth.access_RieszSolver_A().get_testprec();


            for (std::size_t i=0; i<Qa; ++i) {
                FillWithZeros(Lambda, temp[i]);
                (*A[i]).eval(bf, temp[i]);

                // left preconditioning
                for (auto& el : temp[i]) {
                    el.second *= P(el.first)*(-1.);
                }

                THRESH_NoCopy(temp[i], ETA);
                total += temp[i].size();
            }

            this->A_representors.push_back(temp);
            std::cout << "Current total size of A and F is " << total << std::endl;
        }

    public:
        void
        calc_rb_data()
        {
            // F_F_norms
            std::size_t Qf = this->rb_system.Q_f();
            this->rb_system.F_F_representor_norms.engine().resize((int) Qf, (int) Qf);
            for(std::size_t qf1 = 1; qf1 <= Qf; ++qf1) {
                for (std::size_t qf2 = qf1; qf2 <= Qf; ++qf2) {

                    this->rb_system.F_F_representor_norms(qf1,qf2) =
                          dotProduct(this->F_representors[qf1-1],
                                     this->F_representors[qf2-1]);

                    if(qf1 != qf2) {
                        this->rb_system.F_F_representor_norms(qf2,qf1) = this->rb_system.F_F_representor_norms(qf1,qf2);
                    }
                }
            }

            // A_A_norms
            std::size_t Qa      = this->rb_system.Q_a();
            std::size_t n_bf    = this->n_bf();
            this->rb_system.A_A_representor_norms.resize(n_bf);
            for(std::size_t n1 = 0; n1 < n_bf; ++n1) {
                this->rb_system.A_A_representor_norms[n1].clear();

                for (std::size_t n2 = n1; n2 < n_bf; ++n2) {
                    typename RB_Model::FullColMatrixT A_n1_n2(Qa, Qa);
                    for(std::size_t qa1 = 1; qa1 <= Qa; ++qa1) {
                        for(std::size_t qa2 = qa1; qa2 <= Qa; ++qa2) {
                            A_n1_n2(qa1, qa2) =
                            dotProduct(this->A_representors[n1][qa1-1],
                                       this->A_representors[n2][qa2-1]);

                            if(qa1 != qa2){
                                if(n1 == n2){
                                    A_n1_n2(qa2, qa1) = A_n1_n2(qa1, qa2);
                                }
                                else{
                                    A_n1_n2(qa2, qa1) =
                                    dotProduct(this->A_representors[n1][qa2-1],
                                               this->A_representors[n2][qa1-1]);
                                }
                            }
                        }
                    }

                    this->rb_system.A_A_representor_norms[n1].push_back(A_n1_n2);
                }
            }

            // A_F_norms
            this->rb_system.A_F_representor_norms.resize(n_bf);
            for(std::size_t n = 0; n < n_bf; ++n) {
            typename RB_Model::FullColMatrixT A_F(Qa, Qf);
                for(std::size_t qa = 1; qa <= Qa; ++qa) {
                    for(std::size_t qf = 1; qf <= Qf; ++qf) {
                        A_F(qa, qf) =
                        dotProduct(this->A_representors[n][qa-1],
                                   this->F_representors[qf-1]);
                    }
                }
                this->rb_system.A_F_representor_norms[n] = A_F;
            }

            // RB_A
            this->rb_system.RB_A_matrices.resize(Qa);
            for (std::size_t q_a = 0; q_a < Qa; ++q_a) {
                this->rb_system.RB_A_matrices[q_a].resize(n_bf, n_bf);
                for (std::size_t i = 1; i <= n_bf; ++i) {
                    for (std::size_t j = 1; j <= n_bf; ++j) {
                        this->rb_system.RB_A_matrices[q_a](j,i) = this->rb_truth.lhs_u_u(q_a, this->rb_basisfunctions[j-1],
                                                                                        this->rb_basisfunctions[i-1]);
                    }
                }
                std::cout << std::endl << "||------- RB_A (" << q_a << ")  -----------------------------------------||" << std::endl;
                std::cout << this->rb_system.RB_A_matrices[q_a] << std::endl;
            }

            // RB_F
            this->rb_system.RB_F_vectors.resize(Qf);
            for (std::size_t q_f = 0; q_f < Qf; ++q_f) {
                this->rb_system.RB_F_vectors[q_f].resize(n_bf);
                for (std::size_t n = 1; n <= n_bf; ++n) {
                    this->rb_system.RB_F_vectors[q_f](n) = this->rb_truth.rhs_u(q_f, this->rb_basisfunctions[n-1]);
                }
                std::cout << std::endl << "||------- RB_F (" << q_f << ")  -----------------------------------------||" << std::endl;
                std::cout << this->rb_system.RB_F_vectors[q_f] << std::endl;
            }

            // RB_InnerProduct
            this->rb_system.RB_inner_product.resize(n_bf, n_bf);
            for (unsigned int i = 1; i <= n_bf; ++i) {
                for (unsigned int j = i; j <= n_bf; ++j) {
                    this->rb_system.RB_inner_product(i, j) = this->rb_truth.innprod_Y_u_u(this->rb_basisfunctions[i-1],
                                                                                    this->rb_basisfunctions[j-1]);
                    this->rb_system.RB_inner_product(j, i) = this->rb_system.RB_inner_product(i,j);
                }
            }

            std::cout << std::endl << "||------- RB_InnerProduct  ----------------------------------||" << std::endl;
            std::cout << this->rb_system.RB_inner_product << std::endl;
        }
};

        //===============================================================//
        //========= Simple RB System	  ===============================//
        //===============================================================//

        template <typename _T, typename _ParamType, typename LB_System>
        class Simple_RB_System: public RB_System<_T,_ParamType> {

            typedef std::map<_ParamType, _T>   StabConstVec;

        public:
            Simple_RB_System(ThetaStructure<_ParamType>& _thetas_a,
                          ThetaStructure<_ParamType>& _thetas_f,
                          LB_System& _lb_system)
             : RB_System<_T,_ParamType>(_thetas_a, _thetas_f),
               lb_system(_lb_system){}

            _T
            alpha_LB(_ParamType& mu)
    {
    	/*// If it is precalculated, return value
    	for(auto& el: alphas){
    	    bool is_mu = true;
			for(std::size_t i = 0; i < mu.size(); ++i){
				if( (mu[i]-(el.first)[i]) > 1e-4 ){
					is_mu = false;
					break;
				}
			}
			if(is_mu){
				return el.second;
			}
    	}
    	// Else we have to calculate it
    	_T alpha = lb_system.calculate_alpha(mu);
    	alphas.insert(std::make_pair(mu, alpha));
    	return alpha;
    	*/
        return std::min(mu[0], 1.);
    }

    void
    read_alpha(const char* filename){
        ifstream alphaFile(filename);
        while(alphaFile.good()){
             _ParamType mu;
             _T alpha;
             size_t pdim = ParamInfo<_ParamType>::dim;
             for(size_t p = 0; p < pdim; ++p){
            	 alphaFile >> mu[p];
             }
             alphaFile >> alpha;
             alphas.insert(std::make_pair(mu, alpha));
        }
    }

    void
    write_alpha(const char* filename){
        size_t pdim = ParamInfo<_ParamType>::dim;
        ofstream alphaFile(filename);
        if(alphaFile.is_open()){
        	for(auto& a : alphas){
                for(size_t p = 0; p < pdim; ++p){
               	 alphaFile << a.first[p] << " ";
                }
                alphaFile << a.second << std::endl;
        	}
        }
        else{
        	std::cerr << "Couldn't open file " << filename << " for writing." << std::endl;
        	return;
        }
    }
    
    void
    set_alpha(const StabConstVec& _alphas){
        alphas = _alphas;
    }

private:
    StabConstVec	alphas;
    LB_System&		lb_system;
};


//===============================================================//
//========= TYPEDEFS  =======================//
//===============================================================//

//==== General ====//
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef Coefficients<Lexicographical,T,Index2D>						CoeffVector;

//==== Basis 1D & 2D ====//
const DomainType      DomainType_XY                                 = Interval;
const FunctionSide    FctSide_X                                     = Primal;
const Construction    Constr_X                                      = Dijkema;
const Construction    RefConstr_X                                   = Dijkema;
const FunctionSide    FctSide_Y                                     = Orthogonal;
const Construction    Constr_Y                                      = Multi;
const Construction    RefConstr_Y                                   = MultiRefinement;
typedef Basis<T,FctSide_X, DomainType_XY, Constr_X>					Basis_X;
typedef Basis<T,FctSide_Y, DomainType_XY, Constr_Y>					Basis_Y;
typedef Basis_X::RefinementBasis                                    RefBasis_X;
typedef Basis_Y::RefinementBasis                           	        RefBasis_Y;

typedef TensorBasis2D<Adaptive,Basis_X,Basis_Y>                     Basis2D;

//==== Adaptive Operators ====//
typedef AdaptiveIdentityOperator1D<T,FctSide_X, DomainType_XY, Constr_X>                Identity1D_X;
typedef AdaptiveIdentityOperator1D<T,FctSide_Y, DomainType_XY, Constr_Y>                Identity1D_Y;
typedef AdaptiveLaplaceOperator1D<T,FctSide_X, DomainType_XY, Constr_X>                 Laplace1D_X;
typedef AdaptiveLaplaceOperator1D<T,FctSide_Y, DomainType_XY, Constr_Y>                 Laplace1D_Y;
typedef AdaptiveWeightedPDEOperator1D<T,FctSide_X, DomainType_XY, Constr_X>             WeightedLaplace1D_X;
typedef AdaptiveWeightedPDEOperator1D<T,FctSide_X, DomainType_XY, Constr_X>             WeightedIdentity1D_X;

typedef AdaptiveIdentityOperator1D<T,FctSide_X, DomainType_XY, RefConstr_X>             RefIdentity1D_X;
typedef AdaptiveIdentityOperator1D<T,FctSide_Y, DomainType_XY, RefConstr_Y>             RefIdentity1D_Y;
typedef AdaptiveLaplaceOperator1D<T,FctSide_X, DomainType_XY, RefConstr_X>              RefLaplace1D_X;
typedef AdaptiveLaplaceOperator1D<T,FctSide_Y, DomainType_XY, RefConstr_Y>              RefLaplace1D_Y;
typedef AdaptiveWeightedPDEOperator1D<T,FctSide_X, DomainType_XY, RefConstr_X>          RefWeightedLaplace1D_X;
typedef AdaptiveWeightedPDEOperator1D<T,FctSide_X, DomainType_XY, RefConstr_X>          RefWeightedIdentity1D_X;

//==== LocalOperators ====//
typedef LocalOperator1D<Basis_X,Basis_X,RefIdentity1D_X,Identity1D_X>                   LOp_Id1D_X;
typedef LocalOperator1D<Basis_Y,Basis_Y,RefIdentity1D_Y,Identity1D_Y>                   LOp_Id1D_Y;
typedef LocalOperator1D<Basis_X,Basis_X,RefLaplace1D_X,Laplace1D_X>                     LOp_Lapl1D_X;
typedef LocalOperator1D<Basis_Y,Basis_Y,RefLaplace1D_Y,Laplace1D_Y>                     LOp_Lapl1D_Y;
typedef LocalOperator1D<Basis_X,Basis_X,RefWeightedLaplace1D_X,WeightedLaplace1D_X>     LOp_WLapl1D_X;
typedef LocalOperator1D<Basis_X,Basis_X,RefWeightedIdentity1D_X,WeightedIdentity1D_X>   LOp_WId1D_X;

typedef LocalOperator2D<LOp_Id1D_X, LOp_Id1D_Y>                        	                LOp_Id_Id_2D;
typedef LocalOperator2D<LOp_Lapl1D_X, LOp_Id1D_Y>                                       LOp_Lapl_Id_2D;
typedef LocalOperator2D<LOp_Id1D_X, LOp_Lapl1D_Y>                                       LOp_Id_Lapl_2D;
typedef SimpleOp2D<LOp_WLapl1D_X, LOp_Id1D_Y>                                      LOp_WLapl_Id_2D;
typedef SimpleOp2D<LOp_WId1D_X, LOp_Lapl1D_Y>                                      LOp_WId_Lapl_2D;

//==== CompoundOperators ====//
typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 		        Flex_COp_2D;
typedef CompoundLocalOperator<Index2D,LOp_Id_Id_2D,LOp_Lapl_Id_2D,LOp_Id_Lapl_2D>       H1_InnProd_2D;

//==== Preconditioners ====//
typedef H1NormPreconditioner2D<T, Basis2D>                          Prec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//==== RightHandSides ====//
typedef SeparableRHS2D<T,Basis2D>                                   SeparableRhsIntegral2D;
typedef SimpleRHS<T,Index2D,SeparableRhsIntegral2D,
            NoPrec2D>                                         		SeparableRhs;
typedef SimpleRHS<T,Index2D,SeparableRhsIntegral2D,
            Prec2D>                                         		RhsPrec;


//==== RB Stuff ====//
const size_t PDim = 2;

typedef CoeffVector 								 						DataType;
typedef array<T,PDim>	 													ParamType;

typedef AffineLocalOperator<Index2D,AbstractLocalOperator2D<T>,ParamType>	Affine_Op_2D;
typedef AffineRhs<T,Index2D,SeparableRhs,ParamType>							Affine_Rhs_2D;
typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs>							RieszF_Rhs_2D;
typedef FlexibleBilformRhs<Index2D,AbstractLocalOperator2D<T> >				RieszA_Rhs_2D;

typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs>     					Flex_Rhs_2D;

typedef AffineBilformRhs<Index2D, AbstractLocalOperator2D<T>, ParamType>			AffineA_Rhs_2D;
typedef ResidualRhs<Index2D, AffineA_Rhs_2D, Affine_Rhs_2D, ParamType, DataType>	RieszRes_Rhs_2D;


typedef MultiTreeAWGM2<Index2D,Basis2D, Affine_Op_2D,
		Affine_Rhs_2D,Prec2D>                           					MT_AWGM_Truth;
typedef MultiTreeAWGM2<Index2D,Basis2D, H1_InnProd_2D,
		RieszF_Rhs_2D,Prec2D>											    MT_AWGM_Riesz_F;
typedef MultiTreeAWGM2<Index2D,Basis2D, H1_InnProd_2D,
		RieszA_Rhs_2D,Prec2D>											    MT_AWGM_Riesz_A;
typedef MultiTreeAWGM2<Index2D,Basis2D, H1_InnProd_2D,
		RieszRes_Rhs_2D,Prec2D>												MT_AWGM_Riesz_Res;

typedef MT_Truth<DataType,ParamType,MT_AWGM_Truth,
				 MT_AWGM_Riesz_F,MT_AWGM_Riesz_A,H1_InnProd_2D,
				 Flex_COp_2D, Flex_Rhs_2D, MT_AWGM_Riesz_Res>				MTTruthSolver;

typedef LB_Base<ParamType, MTTruthSolver> 									        LB_Type;
typedef Simple_RB_System<T,ParamType, LB_Type>								        RB_Model;
typedef RB_Base_wavest<RB_Model,MTTruthSolver,DataType,ParamType,IndexSet<Index2D>> RB_BaseModel;

#define x1 1./3.
#define x2 2./3.
#define y1 2./5.
#define y2 4./5.
#define nb_stempel 9.

T f_1_x(T x)
{
    if(x <= x1){
        return 1.;
    }
    return 0.;
}

T f_2_x(T x)
{
    if(x >= x1 && x <= x2){
        return 1.;
    }
    return 0.;
}

T f_3_x(T x)
{
    if(x >= x2){
        return 1.;
    }
    return 0.;
}

T f_1_y(T y)
{
    if(y <= y1){
        return 1.;
    }
    return 0.;
}

T f_2_y(T y)
{
    if(y >= y1 && y <= y2){
        return 1.;
    }
    return 0.;
}

T f_3_y(T y)
{
    if(y >= y2){
        return 1.;
    }
    return 0.;
}


T dummy(T, T)
{
    return 0;
}

T no_theta(const std::array<T,PDim>& /*mu*/)
{
	return 1.;
}

T theta_1(const std::array<T,PDim>& mu)
{
	return mu[0];
}

T theta_chi_1(const std::array<T,PDim>& mu)
{
    if(mu[1] <= 1./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_2(const std::array<T,PDim>& mu)
{
    if(mu[1] > 1./nb_stempel && mu[1] <= 2./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_3(const std::array<T,PDim>& mu)
{
    if(mu[1] > 2./nb_stempel && mu[1] <= 3./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_4(const std::array<T,PDim>& mu)
{
    if(mu[1] > 3./nb_stempel && mu[1] <= 4./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_5(const std::array<T,PDim>& mu)
{
    if(mu[1] > 4./nb_stempel && mu[1] <= 5./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_6(const std::array<T,PDim>& mu)
{
    if(mu[1] > 5./nb_stempel && mu[1] <= 6./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_7(const std::array<T,PDim>& mu)
{
    if(mu[1] > 6./nb_stempel && mu[1] <= 7./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_8(const std::array<T,PDim>& mu)
{
    if(mu[1] > 7./nb_stempel && mu[1] <= 8./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_9(const std::array<T,PDim>& mu)
{
    if(mu[1] > 8./nb_stempel){
       return 1;
    }
	return 0;
}

T zero_fct(T /*x*/){
	return 0;
}

T chi_omega_1(T x){
    if(x <= 0.5){
        return 1;
    }
	return 0;
}

T chi_omega_0(T x){
    if(x >= 0.5){
        return 1;
    }
	return 0;
}
