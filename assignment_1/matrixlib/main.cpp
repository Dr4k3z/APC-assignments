#include <iostream>
#include <string>
#include <chrono>
#include <functional>

#include "headers/squarematrix.h"
#include "headers/vectorhelpers.h"
#include "headers/power_iteration.h"
#include "headers/inverse_power_iteration.h"
#include "headers/shift_inverse_power_iteration.h"

void unit_test_max(){
    std::string filenames[] = {"UNIT_TEST/input_10.txt", "UNIT_TEST/input_20.txt", "UNIT_TEST/input_50.txt"};
    eigenvalue::power_iteration pi(10000, 1e-6, BOTH);

    for (auto name : filenames){
        linear_algebra::square_matrix A(name);
        std::vector<double> x0(A.size());
        x0[0] = 1.; //starting point

        double max_expected;

        if (name == "UNIT_TEST/input_10.txt"){
            max_expected = 5.10274;
        }
        else if (name == "UNIT_TEST/input_20.txt"){
            max_expected = 10.1993;
        }
        else if (name == "UNIT_TEST/input_50.txt"){
            max_expected = 25.0613;
        }

        double max_obtained = pi.solve(A, x0);
        std::string result;
        
        if (std::abs(max_obtained - max_expected) < 1e-3)
            result = "CONVERGED";
        else
            result = "NOT CONVERGED";
            
        std::cout << "Maximum modulus eigenvalue: " << max_obtained << " --> " << result << std::endl;
    }
    std::cout << "-----------" << std::endl;
}

void unit_test_min(){
    std::string filenames[] = {"UNIT_TEST/input_10.txt", "UNIT_TEST/input_20.txt", "UNIT_TEST/input_50.txt"};
    eigenvalue::inverse_power_iteration inv_pi(10000, 1e-6, BOTH);
    
    for (auto name : filenames){
        linear_algebra::square_matrix A(name);
        std::vector<double> x0(A.size());
        x0[0] = 1.; //starting point

        double min_expected;

        if (name == "UNIT_TEST/input_10.txt"){
            min_expected = -0.0817798;
        }
        else if (name == "UNIT_TEST/input_20.txt"){
            min_expected = 0.0591506;
        }
        else if (name == "UNIT_TEST/input_50.txt"){
            min_expected = 0.0171826;
        }

        double max_obtained = inv_pi.solve(A, x0);
        std::string result;

        if (std::abs(max_obtained - min_expected) < 1e-3)
            result = "CONVERGED";
        else
            result = "NOT CONVERGED";
    
        std::cout << "Minimum modulus eigenvalue: " << max_obtained << " --> " << result << std::endl;
    }
    std::cout << "-----------" << std::endl;
}

void unit_test_closest(){
    std::string filenames[] = {"UNIT_TEST/input_10.txt", "UNIT_TEST/input_20.txt", "UNIT_TEST/input_50.txt"};
    eigenvalue::shift_inverse_power_iteration s_inv_pi(10000, 1e-6, BOTH);

    for (auto name : filenames){
        linear_algebra::square_matrix A(name);
        std::vector<double> x0(A.size());
        x0[0] = 1.; //starting point

        double closest_expected, mu;

        if (name == "UNIT_TEST/input_10.txt"){
            closest_expected = 0.946734;
            mu = 1;
        }
        else if (name == "UNIT_TEST/input_20.txt"){
            closest_expected = 1.70976;
            mu = 2;
        }
        else if (name == "UNIT_TEST/input_50.txt"){
            closest_expected = 2.04765;
            mu = 2;
        }

        double closest_obtained = s_inv_pi.solve(A, mu, x0);
        std::string result;

        if (std::abs(closest_obtained - closest_expected) < 1e-3)
            result = "CONVERGED";
        else
            result = "NOT CONVERGED";
        std::cout << "Eigenvalue closest to " << mu << ": " << closest_obtained << " --> " << result << std::endl;
    }
    std::cout << "-----------" << std::endl;
}

void time_elapsed(std::function<void()> f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time elapsed: " << std::chrono::duration<double>(end - start).count() << "s" << std::endl;
}

int main() {
    time_elapsed(unit_test_max);
    time_elapsed(unit_test_min);
    time_elapsed(unit_test_closest);
    return 0;
}
