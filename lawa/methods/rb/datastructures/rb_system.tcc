#include <cfloat>
#include <cassert>
#include <limits>
#include <math.h>
#include <iomanip>
#include <sys/stat.h>
#include <lawa/flensforlawa.h>
#include <lawa/settings/enum.h>

namespace lawa {

template<typename T, typename ParamType>
RB_System<T,ParamType>::RB_System(ThetaStructure<ParamType>& _thetas_a,
							 ThetaStructure<ParamType>& _thetas_f)
 : thetas_a(_thetas_a), thetas_f(_thetas_f)
{}

template<typename T, typename ParamType>
std::size_t
RB_System<T,ParamType>::Q_f()
{
	return thetas_f.size();
}

template<typename T, typename ParamType>
std::size_t
RB_System<T,ParamType>::Q_a()
{
	return thetas_a.size();
}

template<typename T, typename ParamType>
flens::DenseVector<flens::Array<T> >
RB_System<T,ParamType>::
get_rb_solution(std::size_t N, ParamType& mu)
{
    using flens::_;

	if(N==0){
		return DenseVectorT(0);
	}

	assert(RB_F_vectors.size() > 0);
	assert(RB_F_vectors[0].length() >= (FLENS_DEFAULT_INDEXTYPE)N);
	assert(RB_A_matrices[0].numRows() >= (FLENS_DEFAULT_INDEXTYPE)N);

	FullColMatrixT A(N, N);

	for (std::size_t i = 0; i < thetas_a.size(); ++i) {
		//A += thetas_a.eval(i,mu) * RB_A_matrices[i](_(1,N), _(1,N));
		flens::blas::axpy(cxxblas::NoTrans, thetas_a.eval(i,mu), RB_A_matrices[i](_(1,N), _(1,N)), A);
	}

	DenseVectorT F(N);
	for (unsigned FLENS_DEFAULT_INDEXTYPE i = 0; i < thetas_f.size(); ++i) {
		//F += thetas_f.eval(i,mu) * RB_F_vectors[i](_(1,N));
		flens::blas::axpy(thetas_f.eval(i,mu), RB_F_vectors[i](_(1,N)), F);
	}

	DenseVectorT u(N);

	FLENS_DEFAULT_INDEXTYPE its;
	switch (rb_params.call) {
		case call_cg:
			its = cg(A, u, F, sqrt(std::numeric_limits<double>::epsilon()), 100*u.length());
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " cg iterations" << std::endl;
			}
			break;
		case call_gmres:
			std::cout << "A = " << A << std::endl;
			std::cout << "F = " << F << std::endl;
			its = gmres(A, u, F, sqrt(std::numeric_limits<double>::epsilon()), 100*u.length());
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " gmres iterations" << std::endl;
                DenseVectorT Au = A*u;
                std::cout << "  A * u = " << Au << std::endl;
                Au -= F;
                std::cout << "  A*u - F = " << std::setprecision(20) << Au << std::endl;
			}
			break;
		case call_cgls:
		    std::cout << "A = " << A << std::endl;
			std::cout << "F = " << F << std::endl;
            std::cout << "tol = " << std::numeric_limits<T>::epsilon() << std::endl;
			its = cgls(A, u, F);
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " cgls iterations" << std::endl;
				DenseVectorT Au = A*u;
                std::cout << "  A * u = " << Au << std::endl;
                Au -= F;
                std::cout << "  A*u - F = " << std::setprecision(20) << Au << std::endl;
			}
			break;	
		default:
			if(rb_params.verbose){
				std::cerr << "RB solution: unrecognized solver call " << std::endl;
			}
			exit(1);
			break;
	}

	return u;
}

template<typename T, typename ParamType>
flens::DenseVector<flens::Array<T> >
RB_System<T,ParamType>::
get_rb_solution(std::vector<std::size_t> indices, ParamType& mu)
{
    using flens::_;

	std::size_t N = indices.size();
	if(N==0){
		return DenseVectorT(0);
	}

	// First assemble all up to largest index...
	std::size_t Nmax = *(std::max_element(indices.begin(), indices.end()))+1;

	assert(RB_F_vectors.size() > 0);
	assert(RB_F_vectors[0].length() >= (FLENS_DEFAULT_INDEXTYPE)Nmax);
	assert(RB_A_matrices[0].numRows() >= (FLENS_DEFAULT_INDEXTYPE)Nmax);

	FullColMatrixT Afull(Nmax, Nmax);

	for (std::size_t i = 0; i < thetas_a.size(); ++i) {
		//A += thetas_a.eval(i,mu) * RB_A_matrices[i](_(1,N), _(1,N));
		flens::blas::axpy(cxxblas::NoTrans, thetas_a.eval(i,mu), RB_A_matrices[i](_(1,Nmax), _(1,Nmax)), Afull);
	}

	DenseVectorT Ffull(Nmax);
	for (std::size_t i = 0; i < thetas_f.size(); ++i) {
		//F += thetas_f.eval(i,mu) * RB_F_vectors[i](_(1,N));
		flens::blas::axpy(thetas_f.eval(i,mu), RB_F_vectors[i](_(1,Nmax)), Ffull);
	}

	// Then throw away unnecessary info
	FullColMatrixT  A(N,N);
	DenseVectorT	F(N);

	std::size_t count_i=1;
	for(auto& i : indices){
		F(count_i) = Ffull(i+1);
		//std::cout << F << std::endl;
		std::size_t count_j=1;
		for(auto& j : indices){
			A(count_i, count_j) = Afull(i+1,j+1);
			count_j++;
		}
		count_i++;
	}

	// Then use the rest
	DenseVectorT u(N);

	FLENS_DEFAULT_INDEXTYPE its;
	switch (rb_params.call) {
		case call_cg:
			its = cg(A, u, F);
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " cg iterations" << std::endl;
			}
			break;
		case call_gmres:
			std::cout << "A = " << A << std::endl;
			std::cout << "F = " << F << std::endl;
			its = gmres(A, u, F);
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " gmres iterations" << std::endl;
				DenseVectorT Au = A*u;
                std::cout << "  A * u = " << Au << std::endl;
                Au -= F;
                std::cout << "  A*u - F = " << std::setprecision(20) << Au << std::endl;
			}
			break;
		case call_cgls:
		    std::cout << "A = " << A << std::endl;
			std::cout << "F = " << F << std::endl;
            std::cout << "tol = " << std::numeric_limits<T>::epsilon() << std::endl;
			its = cgls(A, u, F);
			if(rb_params.verbose){
				std::cout << "RB solution: " << its << " cgls iterations" << std::endl;
				DenseVectorT Au = A*u;
                std::cout << "  A * u = " << Au << std::endl;
                Au -= F;
                std::cout << "  A*u - F = " << std::setprecision(20) << Au << std::endl;
			}
			break;
		default:
			if(rb_params.verbose){
				std::cerr << "RB solution: unrecognized solver call " << std::endl;
			}
			exit(1);
			break;
	}

	return u;
}



template<typename T, typename ParamType>
void
RB_System<T,ParamType>::
get_rb_LGS(std::size_t N, ParamType& mu, FullColMatrixT& A, DenseVectorT& F){
    using flens::_;

    FullColMatrixT A_(N, N);
	for (std::size_t i = 0; i < thetas_a.size(); ++i) {
		//A += thetas_a.eval(i,mu) * RB_A_matrices[i](_(1,N), _(1,N));
		flens::blas::axpy(cxxblas::NoTrans, thetas_a.eval(i,mu), RB_A_matrices[i](_(1,N), _(1,N)), A);
	}

    DenseVectorT F_(N);
	for (unsigned FLENS_DEFAULT_INDEXTYPE i = 0; i < thetas_f.size(); ++i) {
		//F += thetas_f.eval(i,mu) * RB_F_vectors[i](_(1,N));
		flens::blas::axpy(thetas_f.eval(i,mu), RB_F_vectors[i](_(1,N)), F);
	}
	
    flens::blas::copy(F, F_);
    flens::blas::copy(cxxblas::NoTrans, A, A_);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
get_errorbound(const DenseVectorT& u_N, ParamType& mu)
{
    return  residual_dual_norm(u_N, mu) / alpha_LB(mu);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
get_errorbound(std::vector<std::size_t> indices, const DenseVectorT& u_N, ParamType& mu)
{
    return  residual_dual_norm(indices, u_N, mu) / alpha_LB(mu);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
get_errorbound_accuracy(const DenseVectorT& u_N, ParamType& mu,  std::vector<T> eps_f, std::vector<std::vector<T> > eps_a){

	T val = 0;
    for (std::size_t i = 1; i <= thetas_f.size(); ++i) {
        for (std::size_t j = 1; j <= thetas_f.size(); ++j) {
        	//std::cout << "F_F (q=" << i << ",q'=" << j << ") = " << std::max(thetas_f.eval(i-1,mu)*thetas_f.eval(j-1, mu), 0.) << " * " << eps_f << " * " << eps_f << std::endl;
        	val += std::max(thetas_f.eval(i-1,mu)*thetas_f.eval(j-1, mu), 0.) * eps_f[i-1] * eps_f[j-1];
        }
    }

    for(FLENS_DEFAULT_INDEXTYPE i = 1; i <= u_N.length(); ++i){
        for(FLENS_DEFAULT_INDEXTYPE j = 1; j <= u_N.length(); ++j){
            for (std::size_t k = 1; k <= thetas_a.size(); ++k) {
                for (std::size_t l = 1; l <= thetas_a.size(); ++l) {

                	//std::cout << "A_A (n=" << i << ",n'=" << j << ",q=" << k << ",q'=" << l << ") = "
                	//		  << std::max(u_N(i)*u_N(j)*thetas_a.eval(k-1,mu)*thetas_a.eval(l-1,mu),0.) << " * " << eps_a << " * " << eps_a << std::endl;

                	val += std::max(u_N(i)*u_N(j)*thetas_a.eval(k-1,mu)*thetas_a.eval(l-1,mu),0.) * eps_a[i-1][k-1] * eps_a[j-1][l-1];
                }
            }
        }
    }

    for(FLENS_DEFAULT_INDEXTYPE i = 1; i <= u_N.length(); ++i){
        for (std::size_t j = 1; j <= thetas_f.size(); ++j) {
            for (std::size_t k = 1; k <= thetas_a.size(); ++k) {
            	//std::cout << "A_F (n=" << i << ",qf=" << j << ",qa=" << k << ") = "
            	//		  << 2.* std::max(u_N(i)*thetas_f.eval(j-1, mu)*thetas_a.eval(k-1,mu), 0.) << " * " << eps_a << " * " << eps_f << std::endl;

            	val += 2.* std::max(u_N(i)*thetas_f.eval(j-1, mu)*thetas_a.eval(k-1,mu), 0.) * eps_f[j-1] * eps_a[i-1][k-1];
            }
        }
    }

    return std::sqrt(val);
}


template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
residual_dual_norm(const DenseVectorT& u_N, ParamType& mu)
{

    T res_dual_norm = 0;

    FLENS_DEFAULT_INDEXTYPE N = u_N.length();

    std::size_t Qf = thetas_f.size();
    std::size_t Qa = thetas_a.size();

	assert(F_F_representor_norms.numRows() >= (FLENS_DEFAULT_INDEXTYPE)Qf);
	assert(A_F_representor_norms.size() >= (std::size_t)N);
	assert(A_A_representor_norms.size() >= (std::size_t)N);


    DenseVectorT ThetaF(Qf);
    DenseVectorT ThetaA(Qa);
    for (std::size_t i = 1; i <= Qf; ++i) {
        ThetaF(i) = thetas_f.eval(i-1,mu);
    }
    //std::cout << "ThetaF = " << ThetaF << std::endl;
    for (std::size_t i = 1; i <= Qa; ++i) {
        ThetaA(i) = thetas_a.eval(i-1,mu);
    }

    DenseVectorT FF_T = F_F_representor_norms * ThetaF;
    //std::cout << "FF_T = " << FF_T << std::endl;

    res_dual_norm = ThetaF * FF_T;

    //std::cout << "res_dual_norm = " << res_dual_norm << std::endl;


	if(N==0){
		return std::sqrt(res_dual_norm);
	}

    DenseVectorT T_AF_T(N);
    FullColMatrixT T_AA_T(N,N);
    for (FLENS_DEFAULT_INDEXTYPE n1 = 1; n1 <= N; ++n1) {
        DenseVectorT AF_T = A_F_representor_norms[n1-1] * ThetaF;
        T_AF_T(n1) = ThetaA * AF_T;
        for(FLENS_DEFAULT_INDEXTYPE n2 = n1; n2 <= N; ++n2) {
            DenseVectorT AA_T = A_A_representor_norms[n1-1][n2-n1] * ThetaA;
            DenseVectorT T_AA = transpose(A_A_representor_norms[n1-1][n2-n1]) * ThetaA;
            T_AA_T(n1, n2) = ThetaA * AA_T;
            T_AA_T(n2, n1) = ThetaA * T_AA;
        }
    }

    //std::cout << " Residual Dual Norm: size(u) = " << u_RB.length() << ", size(T_AF_T) = " << T_AF_T.length() << std::endl;
    //res_dual_norm += 2 * u_RB * T_AF_T;

    DenseVectorT T_AA_T_u = T_AA_T * u_N;
    res_dual_norm += u_N * T_AA_T_u;
    res_dual_norm += 2. * (u_N * T_AF_T);

    if(res_dual_norm < 0){
      std::cout << "Warning: Residual dual norm negative: " << std::setprecision(20) << res_dual_norm << std::endl;
      res_dual_norm = std::fabs(res_dual_norm);
    }

	std::cout << " res_dual_norm = " << std::sqrt(res_dual_norm) << std::endl;

    return std::sqrt(res_dual_norm);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
residual_dual_norm(std::vector<std::size_t> indices, const DenseVectorT& u_N, ParamType& mu)
{

    T res_dual_norm = 0;

	std::size_t N = indices.size();
	std::size_t Nmax = *(std::max_element(indices.begin(), indices.end()))+1;

    std::size_t Qf = thetas_f.size();
    std::size_t Qa = thetas_a.size();

	assert(F_F_representor_norms.numRows() >= (FLENS_DEFAULT_INDEXTYPE)Qf);
	assert(A_F_representor_norms.size() >= (std::size_t)Nmax);
	assert(A_A_representor_norms.size() >= (std::size_t)Nmax);


    DenseVectorT ThetaF(Qf);
    DenseVectorT ThetaA(Qa);
    for (std::size_t i = 1; i <= Qf; ++i) {
        ThetaF(i) = thetas_f.eval(i-1,mu);
    }
    for (std::size_t i = 1; i <= Qa; ++i) {
        ThetaA(i) = thetas_a.eval(i-1,mu);
    }

    DenseVectorT FF_T = F_F_representor_norms * ThetaF;

    res_dual_norm = ThetaF * FF_T;

	if(N==0){
		return std::sqrt(res_dual_norm);
	}

    DenseVectorT T_AF_T(N);
    FullColMatrixT T_AA_T(N,N);
    std::vector<std::size_t> indices_1based;
    for(auto& i : indices){
    	indices_1based.push_back(i+1);
    }

    std::size_t count_n1=1;
    for (auto& n1 : indices_1based) {
        DenseVectorT AF_T = A_F_representor_norms[n1-1] * ThetaF;
        T_AF_T(count_n1) = ThetaA * AF_T;
        //for(FLENS_DEFAULT_INDEXTYPE n2 = n1; n2 <= N; ++n2) {
        std::size_t count_n2 = 1;
        for(auto& n2 : indices_1based){
        	if(n2 < n1){
        		count_n2++;
        		continue;
        	}

            DenseVectorT AA_T = A_A_representor_norms[n1-1][n2-n1] * ThetaA;
            DenseVectorT T_AA = transpose(A_A_representor_norms[n1-1][n2-n1]) * ThetaA;
            T_AA_T(count_n1, count_n2) = ThetaA * AA_T;
            T_AA_T(count_n2, count_n1) = ThetaA * T_AA;
            count_n2++;
        }
        count_n1++;
    }

    //std::cout << " Residual Dual Norm: size(u) = " << u_RB.length() << ", size(T_AF_T) = " << T_AF_T.length() << std::endl;
    //res_dual_norm += 2 * u_RB * T_AF_T;

    DenseVectorT T_AA_T_u = T_AA_T * u_N;
    res_dual_norm += u_N * T_AA_T_u;
    res_dual_norm += 2. * (u_N * T_AF_T);

    if(res_dual_norm < 0){
      std::cout << "Warning: Residual dual norm negative: " << std::setprecision(20) << res_dual_norm << std::endl;
      res_dual_norm = std::fabs(res_dual_norm);
    }

	//std::cout << " res_dual_norm = " << std::sqrt(res_dual_norm) << std::endl;

    return std::sqrt(res_dual_norm);
}

template<typename T, typename ParamType>
T
RB_System<T,ParamType>::
alpha_LB(ParamType& mu)
{
    T alpha_lb = DBL_MAX;
    for(std::size_t qa = 0; qa < thetas_a.size(); ++qa){
        T reftheta = thetas_a.eval(qa, rb_params.ref_param);
        T mutheta = thetas_a.eval(qa, mu);
        if ((mutheta / reftheta) < alpha_lb) {
            alpha_lb = mutheta / reftheta;
        }
    }

	std::cout << " alpha_LB = " << alpha_lb << std::endl;
    return alpha_lb;
}

template<typename T, typename ParamType>
void
RB_System<T,ParamType>::
write_rb_data(const std::string& directory_name){

	const unsigned FLENS_DEFAULT_INDEXTYPE highprecision = 40;

	// Make a directory to store all the data files
	if(mkdir(directory_name.c_str(), 0777) == -1)
	{
		if(rb_params.verbose){
			  std::cerr << "         [In RB_System::write_rb_data: Directory "
					    << directory_name << " already exists, overwriting contents.]" << std::endl;
		}
	}

	std::size_t n_bf = RB_inner_product.numRows();
	std::string n_bf_filename = directory_name + "/n_bf.txt";
	std::ofstream n_bf_file(n_bf_filename.c_str());
	n_bf_file << n_bf << std::endl;
	n_bf_file.close();

	// Write RB_A_matrices
	for(std::size_t i = 0; i < Q_a(); ++i){
		std::stringstream filename;
		filename << directory_name << "/RB_A_" << i+1 << ".txt";
		std::ofstream file(filename.str().c_str());
		file.precision(highprecision);
		file << std::scientific << RB_A_matrices[i] << std::endl;
		file.close();
	}

	// Write RB_F_vectors
	for(std::size_t i = 0; i < Q_f(); ++i){
		std::stringstream filename;
		filename << directory_name << "/RB_F_" << i+1 << ".txt";
		std::ofstream file(filename.str().c_str());
		file.precision(highprecision);
		file << std::scientific << RB_F_vectors[i] << std::endl;
		file.close();
	}

	// Write RB_inner_product
	std::stringstream filename;
	filename << directory_name << "/RB_inner_product.txt";
	std::ofstream file(filename.str().c_str());
	file.precision(highprecision);
	file << std::scientific << RB_inner_product << std::endl;
	file.close();

	// Write F_F_representor_norms
	std::stringstream repr_F_filename;
	repr_F_filename << directory_name << "/F_F_representor_norms.txt";
	std::ofstream repr_F_file(repr_F_filename.str().c_str());
	repr_F_file.precision(highprecision);
	repr_F_file << std::scientific << F_F_representor_norms << std::endl;
	repr_F_file.close();

	// Write A_F_representor_norms
	std::stringstream repr_A_F_filename;
	repr_A_F_filename << directory_name << "/A_F_representor_norms.txt";
	std::ofstream repr_A_F_file(repr_A_F_filename.str().c_str());
	repr_A_F_file.precision(highprecision);
	for(std::size_t i = 0; i < n_bf; ++i){
		repr_A_F_file << std::scientific << A_F_representor_norms[i] << std::endl;
	}
	repr_A_F_file.close();

 	// Write A_A_representor_norms
	std::stringstream repr_A_A_filename;
	repr_A_A_filename << directory_name << "/A_A_representor_norms.txt";
	std::ofstream repr_A_A_file(repr_A_A_filename.str().c_str());
	repr_A_A_file.precision(highprecision);
	for(std::size_t i = 0; i < n_bf; ++i){
		for(std::size_t j = i; j < n_bf; ++j){
			repr_A_A_file << std::scientific << A_A_representor_norms[i][j-i] << std::endl;
		}
	}
	repr_A_A_file.close();

}

template<typename T, typename ParamType>
void
RB_System<T,ParamType>::
read_rb_data(const std::string& directory_name, FLENS_DEFAULT_INDEXTYPE nb)
{
    using flens::_;

	// Read Nr of BasisFunctions
	unsigned FLENS_DEFAULT_INDEXTYPE n_bf;
	std::string n_bf_filename = directory_name + "/n_bf.txt";

	std::ifstream n_bf_file(n_bf_filename.c_str());
	if(n_bf_file.is_open()){
		n_bf_file >> n_bf;
	    n_bf_file.close();

	    if(rb_params.verbose){
		    std::cout << "Number of Basis Functions: " << n_bf << std::endl;
	    }
	}
	else{
	    std::cerr << "Unable to read number of basis functions: " << strerror(errno) << std::endl;
	    exit(1);
	}

	// If N < 0 or > n_bf, set it to n_bf
	nb = (nb < 0) ? n_bf : nb;
	nb = (nb > (FLENS_DEFAULT_INDEXTYPE)n_bf) ? n_bf : nb;

	// Read RB_A_matrices
	RB_A_matrices.clear();
	for(std::size_t i = 0; i < Q_a(); ++i){
	    std::stringstream filename;
	    filename << directory_name << "/RB_A_" << i+1 << ".txt";
	    std::ifstream file(filename.str().c_str());
	    if(file.is_open()){
	    	FullColMatrixT RB_A(n_bf, n_bf);
	    	for(unsigned FLENS_DEFAULT_INDEXTYPE n1 = 1; n1 <= n_bf; n1++){
	    		for(unsigned FLENS_DEFAULT_INDEXTYPE n2 = 1; n2 <= n_bf; n2++){
	    			file >> RB_A(n1, n2);
	    		}
	    	}
	    	file.close();
	    	RB_A_matrices.push_back(RB_A(_(1,nb), _(1,nb)));
		    if(rb_params.verbose){
		    	std::cout << " Read " << filename.str() << std::endl;
		    }
	    }
	    else{
	    	std::cerr << "Unable to open file " << filename.str() << " for reading: "<< strerror(errno) << std::endl;
	    	exit(1);
	    }
	}

	// Read RB_F_vectors
	RB_F_vectors.clear();
	for(unsigned FLENS_DEFAULT_INDEXTYPE i = 0; i < Q_f(); ++i){
		std::stringstream filename;
	    filename << directory_name << "/RB_F_" << i+1 << ".txt";
	    std::ifstream file(filename.str().c_str());
	    if(file.is_open()){
	    	DenseVectorT RB_F(n_bf);
	    	for(unsigned FLENS_DEFAULT_INDEXTYPE n = 1; n <= n_bf; n++){
	    		file >> RB_F(n);
	    	}
	    	file.close();
	    	RB_F_vectors.push_back(RB_F(_(1,nb)));
		    if(rb_params.verbose){
		    	std::cout << " Read " << filename.str() << std::endl;
		    }
	    }
	    else{
	    	std::cerr << "Unable to open file " << filename.str() << " for reading: " << strerror(errno) << std::endl;
	    	exit(1);
	    }
	 }

	 // Read RB_inner_product
	 RB_inner_product.engine().resize((FLENS_DEFAULT_INDEXTYPE)n_bf, (FLENS_DEFAULT_INDEXTYPE)n_bf);
	 std::stringstream inner_product_filename;
	 inner_product_filename << directory_name << "/RB_inner_product.txt";
	 std::ifstream inner_product_file(inner_product_filename.str().c_str());
	 if(inner_product_file.is_open()){
		 for(unsigned FLENS_DEFAULT_INDEXTYPE n1 = 1; n1 <= n_bf; n1++){
			 for(unsigned FLENS_DEFAULT_INDEXTYPE n2 = 1; n2 <= n_bf; n2++){
				 inner_product_file >> RB_inner_product(n1, n2);
			 }
		 }
		 if(nb < (FLENS_DEFAULT_INDEXTYPE)n_bf){
			 FullColMatrixT tmp(RB_inner_product);
			 RB_inner_product.engine().resize(nb, nb);
			 RB_inner_product = tmp(_(1,nb), _(1,nb));
		 }
		 inner_product_file.close();
		 if(rb_params.verbose){
			 std::cout << " Read " << inner_product_filename.str() << std::endl;
		 }
	 }
	 else{
		 std::cerr << "Unable to open file " << inner_product_filename.str() << " for reading: " << strerror(errno) << std::endl;
		 exit(1);
	 }

	 // Read F_F_representor norms
	 F_F_representor_norms.engine().resize((FLENS_DEFAULT_INDEXTYPE)Q_f(), (FLENS_DEFAULT_INDEXTYPE)Q_f());
	 std::stringstream F_F_filename;
	 F_F_filename << directory_name << "/F_F_representor_norms.txt";
	 std::ifstream F_F_file(F_F_filename.str().c_str());
	 if(F_F_file.is_open()){
		 for(unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= Q_f(); ++i){
			 for(unsigned FLENS_DEFAULT_INDEXTYPE j = 1; j <= Q_f(); ++j){
				 F_F_file >> F_F_representor_norms(i,j);
			 }
		 }
		 F_F_file.close();
		 if(rb_params.verbose){
			 std::cout << " Read " << F_F_filename.str() << std::endl;
		 }
	 }
	 else{
		 std::cerr << "Unable to open file " << F_F_filename.str() << " for reading: " << strerror(errno) << std::endl;
		 exit(1);
	 }

	  /*// Read output_output_representor norms
	  output_output_representor_norms.engine().resize((FLENS_DEFAULT_INDEXTYPE)Q_output(), (FLENS_DEFAULT_INDEXTYPE)Q_output());
	  std::stringstream output_output_filename;
	  output_output_filename << directory_name << "/output_output_representor_norms.dat";
	  std::ifstream output_output_file(output_output_filename.str().c_str());
	  if(output_output_file.is_open()){
	    for(unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= Q_output(); ++i){
	      for(unsigned FLENS_DEFAULT_INDEXTYPE j = 1; j <= Q_output(); ++j){
	        output_output_file >> output_output_representor_norms(i,j);
	      }
	    }
	    output_output_file.close();
	    std::cout << " Read " << output_output_filename.str() << std::endl;
	  }
	  else{
	    std::cerr << "Unable to open file " << output_output_filename.str() << " for reading!" << std::endl;
	    exit(1);
	  }*/

	 // Read A_F_representor norms
	 A_F_representor_norms.clear();
	 std::stringstream A_F_filename;
	 A_F_filename << directory_name << "/A_F_representor_norms.txt";
	 std::ifstream A_F_file(A_F_filename.str().c_str());
	 if(A_F_file.is_open()){
		 for(FLENS_DEFAULT_INDEXTYPE n = 1; n <= nb; ++n){
			 FullColMatrixT A_F_matrix(Q_a(), Q_f());
			 for(unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= Q_a(); ++i){
				 for(unsigned FLENS_DEFAULT_INDEXTYPE j = 1; j <= Q_f(); ++j){
					 A_F_file >> A_F_matrix(i,j);
				 }
			 }
			 A_F_representor_norms.push_back(A_F_matrix);
		 }
		 A_F_file.close();
		 if(rb_params.verbose){
			 std::cout << " Read " << A_F_filename.str() << std::endl;
		 }
	 }
	 else{
		 std::cerr << "Unable to open file " << A_F_filename.str() << " for reading: " << strerror(errno) << std::endl;
		 exit(1);
	 }

	 // Read A_A_representor norms
	 A_A_representor_norms.clear();
	 std::stringstream A_A_filename;
	 A_A_filename << directory_name << "/A_A_representor_norms.txt";
	 std::ifstream A_A_file(A_A_filename.str().c_str());
	 if(A_A_file.is_open()){
		 for(FLENS_DEFAULT_INDEXTYPE n1 = 1; n1 <= nb; ++n1){
			 std::vector<FullColMatrixT> A_A_n_vec;
			 for(unsigned FLENS_DEFAULT_INDEXTYPE n2 = n1; n2 <= n_bf; ++n2){	// We have to read every information...
				 FullColMatrixT A_A_matrix(Q_a(), Q_a());
				 for(unsigned FLENS_DEFAULT_INDEXTYPE i = 1; i <= Q_a(); ++i){
					 for(unsigned FLENS_DEFAULT_INDEXTYPE j = 1; j <= Q_a(); ++j){
						 A_A_file >> A_A_matrix(i,j);
					 }
				 }
				 if((FLENS_DEFAULT_INDEXTYPE)n2 <= nb){				// .. but safe only for n <= nb
					 A_A_n_vec.push_back(A_A_matrix);
				 }
			 }
			 A_A_representor_norms.push_back(A_A_n_vec);
		 }
		 A_A_file.close();
		 if(rb_params.verbose){
			 std::cout << " Read " << A_A_filename.str() << std::endl;
		 }
	 }
	 else{
		 std::cerr << "Unable to open file " << A_A_filename.str() << " for reading: " << strerror(errno) << std::endl;
		 exit(1);
	 }
}


template<typename T, typename ParamType>
void
RB_System<T,ParamType>::
remove_basisfunction(std::size_t nb)
{
    using flens::_;
	FLENS_DEFAULT_INDEXTYPE n_bf = RB_A_matrices[0].numRows();

	if(nb == 1){
		// F_vectors
		for(std::size_t qf = 0; qf < Q_f(); ++qf){
			RB_F_vectors[qf] = RB_F_vectors[qf](_(RB_F_vectors[qf].firstIndex()+1,RB_F_vectors[qf].lastIndex()));
		}
		// A_matrices
		for(std::size_t qa = 0; qa < Q_a(); ++qa){
			RB_A_matrices[qa] = RB_A_matrices[qa](_(RB_A_matrices[qa].firstRow()+1,RB_A_matrices[qa].lastRow()),
												  _(RB_A_matrices[qa].firstCol()+1,RB_A_matrices[qa].lastCol()));
		}
		// Inner product matrix
		RB_inner_product = RB_inner_product(_(RB_inner_product.firstRow()+1,RB_inner_product.lastRow()),
				   	   	   	   	   	   	    _(RB_inner_product.firstCol()+1,RB_inner_product.lastCol()));

	}
	else{
		if((FLENS_DEFAULT_INDEXTYPE)nb == n_bf){
			// F_vectors
			for(std::size_t qf = 0; qf < Q_f(); ++qf){
				RB_F_vectors[qf] = RB_F_vectors[qf](_(RB_F_vectors[qf].firstIndex(),RB_F_vectors[qf].lastIndex()-1));
			}
			// A_matrices
			for(std::size_t qa = 0; qa < Q_a(); ++qa){

				RB_A_matrices[qa] = RB_A_matrices[qa](_(RB_A_matrices[qa].firstRow(),RB_A_matrices[qa].lastRow()-1),
													  _(RB_A_matrices[qa].firstCol(),RB_A_matrices[qa].lastCol()-1));
			}
			// Inner product matrix
			RB_inner_product = RB_inner_product(_(RB_inner_product.firstRow(),RB_inner_product.lastRow()-1),
					   	   	   	   	   	   	    _(RB_inner_product.firstCol(),RB_inner_product.lastCol()-1));
		}
		else{
			// F_vectors
			for(std::size_t qf = 0; qf < Q_f(); ++qf){
				DenseVectorT tmp = RB_F_vectors[qf];
				RB_F_vectors[qf].engine().resize(n_bf-1);
				RB_F_vectors[qf](_(RB_F_vectors[qf].firstIndex(), nb-1)) = tmp(_(tmp.firstIndex(), nb-1));
				RB_F_vectors[qf](_(nb,RB_F_vectors[qf].lastIndex())) = tmp(_(nb+1, tmp.lastIndex()));

			}
			// A_matrices
			for(std::size_t qa = 0; qa < Q_a(); ++qa){
				FullColMatrixT tmp = RB_A_matrices[qa];
				RB_A_matrices[qa].engine().resize(n_bf-1, n_bf-1);
				RB_A_matrices[qa](_(RB_A_matrices[qa].firstRow(),nb-1),_(RB_A_matrices[qa].firstCol(),nb-1))
					= tmp(_(tmp.firstRow(), nb-1), _(tmp.firstCol(), nb-1));
				RB_A_matrices[qa](_(nb, RB_A_matrices[qa].lastRow()),_(RB_A_matrices[qa].firstCol(),nb-1))
					= tmp(_(nb+1, tmp.lastRow()), _(tmp.firstCol(), nb-1));
				RB_A_matrices[qa](_(RB_A_matrices[qa].firstRow(),nb-1),_(nb, RB_A_matrices[qa].lastCol()))
					= tmp(_(tmp.firstRow(), nb-1), _(nb+1, tmp.lastCol()));
				RB_A_matrices[qa](_(nb, RB_A_matrices[qa].lastRow()),_(nb, RB_A_matrices[qa].lastCol()))
					= tmp(_(nb+1, tmp.lastRow()), _(nb+1, tmp.lastCol()));

			}
			// Inner product matrix
			FullColMatrixT tmp = RB_inner_product;
			RB_inner_product.engine().resize(n_bf-1, n_bf-1);
			RB_inner_product(_(RB_inner_product.firstRow(),nb-1), _(RB_inner_product.firstCol(),nb-1))
				= tmp(_(tmp.firstRow(), nb-1), _(tmp.firstCol(), nb-1));
			RB_inner_product(_(nb, RB_inner_product.lastRow()),_(RB_inner_product.firstCol(),nb-1))
				= tmp(_(nb+1, tmp.lastRow()), _(tmp.firstCol(), nb-1));
			RB_inner_product(_(RB_inner_product.firstRow(),nb-1),_(nb, RB_inner_product.lastCol()))
				= tmp(_(tmp.firstRow(), nb-1), _(nb+1, tmp.lastCol()));
			RB_inner_product(_(nb, RB_inner_product.lastRow()),_(nb, RB_inner_product.lastCol()))
				= tmp(_(nb+1, tmp.lastRow()), _(nb+1, tmp.lastCol()));
		}
	}

	// Representor norms

	std::cout << "A_F_representor_norms before: " << std::endl;
	for(auto& el : A_F_representor_norms){
		std::cout << el << std::endl;
	}

	A_F_representor_norms.erase(A_F_representor_norms.begin()+nb-1);
	std::cout << "A_F_representor_norms after: " << std::endl;
	for(auto& el : A_F_representor_norms){
		std::cout << el << std::endl;
	}
	std::cout << std::endl;

	std::cout << "A_A_representor_norms before: " << std::endl;
	for(auto& el : A_A_representor_norms){
		std::cout << "  (next inner loop) " << std::endl;
		for(auto& el2 : el){
			std::cout << el2 << std::endl;
		}
	}
	for(std::size_t  i = 0; i < nb-1; ++i){
		A_A_representor_norms[i].erase(A_A_representor_norms[i].begin()+(nb-i-1));
	}
	A_A_representor_norms.erase(A_A_representor_norms.begin()+nb-1);

	std::cout << "A_A_representor_norms after: " << std::endl;
	for(auto& el : A_A_representor_norms){
		std::cout << "  (next inner loop) " << std::endl;
		for(auto& el2 : el){
			std::cout << el2 << std::endl;
		}
	}
}

} // namespace lawa
