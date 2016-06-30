#include <iostream>
#include <fstream>

namespace lawa {

template<typename T, typename Basis1D>
void 
print_U(const flens::DenseVector<flens::Array<T> >& u, const Basis1D& basis, const FLENS_DEFAULT_INDEXTYPE J, 
        const char* filename, const double deltaX)
{
    std::ofstream file(filename);
    for(double x = 0; x <= 1.; x += deltaX){
        file << x << " " << evaluate(basis,J, u, x, 0) << std::endl; 
    }
    file.close();
}
        
template<typename T, typename Basis2D>
void
print_U(const flens::DenseVector<flens::Array<T> >& u, const Basis2D& basis, const FLENS_DEFAULT_INDEXTYPE J_x, const FLENS_DEFAULT_INDEXTYPE J_y, 
                   const char* filename, const double deltaX, const double deltaY)
{
    std::ofstream file(filename);
    for(double x = 0.; x <= 1.; x += deltaX){
        for(double y = 0; y <= 1.; y += deltaY){
            file << x << " " << y << " " << evaluate(basis, J_x, J_y, u, x, y, 0, 0) << std::endl;
        }
    }
    file.close();
}
        
template<typename T, typename Basis1D>
void 
print_U(const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U, const Basis1D& basis, 
        const FLENS_DEFAULT_INDEXTYPE J, const char* filename, const T timestep, const FLENS_DEFAULT_INDEXTYPE K, const double deltaX)
{
    std::ofstream file(filename);
    for (FLENS_DEFAULT_INDEXTYPE k = 0; k <= K; k++) {
        for(double x = 0; x <= 1.; x += deltaX){
            file << k*timestep << " " << x << " " << evaluate(basis,J, U(U.rows(), k), x, 0) << std::endl; 
        }        
    }

    file.close();
}

} // namespace lawa
