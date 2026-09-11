#include "statistics.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool check(bool condition, const std::string& message)
{
    if (!condition)
    {
        std::cerr << "FAILED: " << message << '\n';
        return false;
    }
    return true;
}

bool nearly_equal(double lhs, double rhs, double tolerance = 1e-10)
{
    return std::abs(lhs - rhs) <= tolerance;
}

} // namespace

int main()
{
    bool passed = true;
    // ---------------------------------------------------------------------
    // Mean
    // ---------------------------------------------------------------------
    std::vector<double> values{1.0, 2.0, 3.0, 4.0, 5.0};
    double mean = mean_value(values);
    passed &= check(nearly_equal(mean, 3.0), "Mean of {1,2,3,4,5} should equal 3");

    // ---------------------------------------------------------------------
    // Population variance
    // ---------------------------------------------------------------------
    double data_variance = variance(values, mean);
    passed &= check(nearly_equal(data_variance, 2.0), "Population variance of {1,2,3,4,5} should equal 2");

    // ---------------------------------------------------------------------
    // Filtering
    //
    // filter = 0.4 removes the first 40% of five values:
    // remaining values = {3,4,5}
    // ---------------------------------------------------------------------
    double filtered_mean = mean_value(values, 0.4);
    passed &= check(nearly_equal(filtered_mean, 4.0), "Filtered mean of {3,4,5} should equal 4");
    double filtered_variance = variance(values, filtered_mean, 0.4);
    passed &= check(nearly_equal(filtered_variance, 2.0 / 3.0), "Filtered population variance should equal 2/3");

    // ---------------------------------------------------------------------
    // Power-law fitting
    //
    // y = 3 * x^2
    // Expected fitted parameters:
    // a = 3
    // b = 2
    // ---------------------------------------------------------------------
    std::vector<double> x_values{1.0, 2.0, 4.0, 8.0};
    std::vector<double> y_values{3.0, 12.0, 48.0, 192.0};
    std::vector<double> parameters = power_fitting(x_values, y_values);
    passed &= check(parameters.size() == 2, "Power-law fitting should return two parameters");

    if (parameters.size() == 2)
    {
        passed &= check(nearly_equal(parameters[0], 3.0), "Power-law coefficient a should equal 3");
        passed &= check(nearly_equal(parameters[1], 2.0), "Power-law exponent b should equal 2");
    }

    if (!passed)
    {
        return 1;
    }
    std::cout << "Statistics unit tests passed.\n";
    return 0;
}
