#include "../headers/inverse_power_iteration.h"

std::vector<double> luSolve(const linear_algebra::square_matrix& U, const linear_algebra::square_matrix& L, const std::vector<double> x){
    std::vector<double> z = linear_algebra::backsolve(U,x);
    return linear_algebra::forwardsolve(L,z);
}

namespace eigenvalue {

    using linear_algebra::operator*;
    using linear_algebra::operator-;

    double inverse_power_iteration::solve(const linear_algebra::square_matrix& A, const std::vector<double>& x0) const{
        /*
            Instead of inverting the matrix, we carry out the following computations.
            Whenever we'd need to compute A^-1*x, we instead solve the system Ay = x, 
            using LU decomposition.
        */
        const float tolerance = 10^-6;
        const unsigned int Tmax = 10000;
        double increment, residual = tolerance+1;

        unsigned int t = 0;
        bool converged = false;

        linear_algebra::square_matrix U,L;
        linear_algebra::lu(A,L,U);
        std::vector<double> x = x0;

        auto Ainvx = luSolve(U,L,x);
        double nu = linear_algebra::scalar(x, Ainvx);

        double nu_new = nu;
        while (!converged && t < Tmax){
            std::vector<double> z = luSolve(U,L,x);
            x = (1/linear_algebra::norm(z))*z;

            auto tmp = luSolve(U,L,x); // to speed things up
            
            nu_new = linear_algebra::scalar(x, tmp);
            residual = linear_algebra::norm(tmp-nu*x);
            increment = abs(nu_new-nu)/nu_new;
            nu = nu_new;   
            if (residual < tolerance && increment < tolerance){
                converged = true;
            }
            t++;
        }

        return 1/nu;
    }

    bool inverse_power_iteration::converged(const double& residual, const double& increment) const
    {
        bool conv;

        switch(termination) {
            case(RESIDUAL):
                conv = residual < tolerance;
                break;
            case(INCREMENT):
                conv = increment < tolerance;
                break;
            case(BOTH):
                conv = residual < tolerance && increment < tolerance;
                break;
            default:
                conv = false;
        }
        return conv;
    }

} // eigenvalue