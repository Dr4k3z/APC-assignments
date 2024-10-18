#ifndef SQUARE_MATRIX_H
#define SQUARE_MATRIX_H

#include <vector>
#include <iostream>
#include <utility>

namespace linear_algebra
{
    class square_matrix
    {
    protected:
        std::vector<double> elem;
        std::size_t sz;

        std::size_t sub2ind(const std::size_t& i, const std::size_t& j) const;
    public:
        square_matrix() = default;
        square_matrix(const std::size_t& s): elem(s*s),  sz(s) {};
        square_matrix(const std::string& filename); //builds a matrix reading data from a text file
        std::size_t size() const { return sz; }

        const double& operator() (const std::size_t& i, const std::size_t& j) const;
        double& operator()(const std::size_t& i, const std::size_t& j);

        square_matrix operator-(const square_matrix& B) const{
            square_matrix C(sz);
            for (std::size_t i = 0; i < sz; i++){
                for (std::size_t j = 0; j < sz; j++){
                    C(i,j) = (*this)(i,j) - B(i,j);
                }
            }
            return C;
        }
        square_matrix operator*(const double scalar) const{
            square_matrix C(sz);
            for (std::size_t i = 0; i < sz; i++){
                for (std::size_t j = 0; j < sz; j++){
                    C(i,j) = (*this)(i,j)*scalar;
                }
            }
            return C;
        }
        std::vector<double> operator*(const std::vector<double>& x) const;
        
        void print(){
            for (std::size_t i = 0; i < sz; i++){
                for (std::size_t j = 0; j < sz; j++){
                    std::cout << elem[sub2ind(i,j)] << " ";
                }
                std::cout << std::endl;
            }
        }
    };

    class identity_matrix : public square_matrix{
    public:
        identity_matrix(const std::size_t& n) : square_matrix(n){
            for (std::size_t i = 0; i < n; i++){
                elem[sub2ind(i,i)] = 1;
            }
        }
    };
}


#endif
