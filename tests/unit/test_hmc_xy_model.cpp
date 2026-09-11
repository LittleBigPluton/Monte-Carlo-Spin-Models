#include "hmc_xy_model.hpp"

#include <cmath>
#include <iostream>
#include <random>
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

class TestHMCXYModel : public HMCXYModel {
public:
    void configure(int length)
    {
        lattice_length = length;
        grid_points = length * length;
        nearest_neighbour();
    }
};

} // namespace

int main()
{
    bool passed = true;
    TestHMCXYModel model;
    model.configure(3);

    // ---------------------------------------------------------------------
    // Fully aligned XY configuration
    //
    // theta_i = 0 for every spin.
    //
    // Every bond contributes -cos(0) = -1.
    // A square periodic lattice contains 2N counted bonds.
    // For N = 9:
    // E = -18
    // |M| = 9
    // ---------------------------------------------------------------------

    std::vector<double> aligned_spins(static_cast<std::size_t>(model.grid_points), 0.0);
    passed &= check(nearly_equal(model.total_energy(aligned_spins), -18.0), "Aligned 3x3 XY configuration should have total energy -18");
    passed &= check(nearly_equal(model.total_magnetization(aligned_spins), 9.0), "Aligned 3x3 XY configuration should have magnetization 9");

    // The derivative of the energy vanishes for a uniform configuration.
    for (int site = 0; site < model.grid_points; ++site)
    {
        passed &= check(nearly_equal(model.gradient(aligned_spins, site), 0.0), "Gradient should vanish for an aligned XY configuration");
    }

    // ---------------------------------------------------------------------
    // Kinetic energy
    //
    // K = 1/2 * sum(p_i^2)
    // ---------------------------------------------------------------------
    std::vector<double> momenta{1.0, -1.0, 2.0, -2.0, 0.0, 0.5, -0.5, 1.5, -1.5};
    double expected_kinetic_energy = 0.0;
    for (double momentum : momenta)
    {
        expected_kinetic_energy += 0.5 * momentum * momentum;
    }

    passed &= check(nearly_equal(model.kinetic_energy(momenta), expected_kinetic_energy), "Kinetic-energy calculation is incorrect");

    // ---------------------------------------------------------------------
    // Hamiltonian
    //
    // H = K + beta * V
    // ---------------------------------------------------------------------
    const double beta = 0.5;
    double expected_hamiltonian =expected_kinetic_energy + beta * model.total_energy(aligned_spins);
    passed &= check(nearly_equal(model.hamiltonian(aligned_spins,momenta,beta), expected_hamiltonian), "Hamiltonian should equal kinetic + beta * potential energy");

    // ---------------------------------------------------------------------
    // Leapfrog equilibrium test
    //
    // Uniform spins + zero momenta gives:
    // gradient = 0
    //
    // Therefore leapfrog should leave the complete state unchanged.
    // ---------------------------------------------------------------------
    std::vector<double> leapfrog_spins = aligned_spins;
    std::vector<double> zero_momenta(static_cast<std::size_t>(model.grid_points),0.0);
    model.leapfrog(zero_momenta, leapfrog_spins, 10, 1.0, 0.05);
    for (int site = 0; site < model.grid_points; ++site)
    {
        passed &= check(nearly_equal(leapfrog_spins[static_cast<std::size_t>(site)], 0.0), "Aligned zero-momentum state should remain unchanged under leapfrog");
        passed &= check(nearly_equal(zero_momenta[static_cast<std::size_t>(site)], 0.0), "Momentum should remain zero for equilibrium configuration");
    }

    // ---------------------------------------------------------------------
    // Random spin initialization
    // ---------------------------------------------------------------------
    std::mt19937 generator(42);
    model.initialize_xy_spins(generator);
    passed &= check(model.xy_spins.size() == static_cast<std::size_t>(model.grid_points), "XY initializer should create one angle per lattice site");
    const double two_pi = 2.0 * std::acos(-1.0);

    for (double angle : model.xy_spins)
    {
        passed &= check(angle >= 0.0 && angle <= two_pi, "Initialized XY angles should lie in [0, 2*pi]");
    }

    if (!passed)
    {
        return 1;
    }
    std::cout << "HMC XY-model unit tests passed.\n";
    return 0;
}
