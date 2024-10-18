#include "../headers/shift_inverse_power_iteration.h"
#include "../headers/inverse_power_iteration.h"
#include "../headers/power_iteration.h"

namespace eigenvalue {

    using linear_algebra::operator*;
    using linear_algebra::operator-;

    double shift_inverse_power_iteration::solve(const linear_algebra::square_matrix& A, const double& mu, const std::vector<double>& x0) const{
        const size_t n = A.size();
        const linear_algebra::identity_matrix id(n);
        
        double lambda_1 = eigenvalue::inverse_power_iteration::solve(A-id*mu,x0);

        return lambda_1+mu;
    }

} // eigenvalue