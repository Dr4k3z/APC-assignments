#include "../headers/power_iteration.h"

namespace eigenvalue {

    using linear_algebra::operator*;
    using linear_algebra::operator-;

    double power_iteration::solve(const linear_algebra::square_matrix& A, const std::vector<double>& x0) const{
        const float tolerance = 10^-6;
        const unsigned int Tmax = 10000;
        double increment, residual = tolerance+1;

        unsigned int t = 0;
        bool converged = false;

        std::vector<double> x = x0;
        double nu = linear_algebra::scalar(x, A*x);
        double nu_new = nu;
        while (!converged && t < Tmax){
            std::vector<double> z = A*x;
            x = (1/linear_algebra::norm(z))*z;

            auto tmp = A*x; // to speed things up
            
            nu_new = linear_algebra::scalar(x, tmp);
            residual = linear_algebra::norm(tmp-nu*x);
            increment = abs(nu_new-nu)/nu_new;
            nu = nu_new;   
            if (residual < tolerance && increment < tolerance){
                converged = true;
            }
            t++;
        }
         
        return nu;
    }

    bool power_iteration::converged(const double& residual, const double& increment) const
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