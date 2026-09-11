#include "ising_model.hpp"
#include "statistics.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    std::mt19937 generator(42);

    // Test the lattice and Ising-model connection.
    Ising model;
    model.initialize();  // Enter a small lattice size, for example 4.
    model.discrete_spin_generator(generator, false, "up");
    model.Temp = 2.3;
    model.total_energy();
    model.total_magnetization();

    std::cout << "Total energy: " << model.Teng << '\n';
    std::cout << "Total magnetization: " << model.Tmag << '\n';

    // Test statistics functions.
    std::vector<double> values{1.0,2.0,3.0,4.0,5.0};
    double average = mean_value(values, 0.0);
    double data_variance = variance(values, average, 0.0);
    std::cout << "Mean: " << average << '\n';
    std::cout << "Variance: " << data_variance << '\n';

    // Simple validation checks.
    if (std::abs(average - 3.0) > 1e-12)
    {
        std::cerr << "Mean test failed.\n";
        return 1;
    }

    if (std::abs(data_variance - 2.0) > 1e-12)
    {
        std::cerr << "Variance test failed.\n";
        return 1;
    }

    std::cout << "All tests completed successfully.\n";

    return 0;
}
