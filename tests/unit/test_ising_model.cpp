#include "ising_model.hpp"

#include <cmath>
#include <iostream>
#include <random>
#include <string>

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

bool nearly_equal(double lhs, double rhs, double tolerance = 1e-12)
{
    return std::abs(lhs - rhs) <= tolerance;
}

class TestIsing : public Ising {
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
    std::mt19937 generator(42);
    TestIsing model;
    model.configure(3);

    // ---------------------------------------------------------------------
    // Cold all-up configuration
    // ---------------------------------------------------------------------
    model.discrete_spin_generator(generator, false, "up");
    passed &= check(model.spins.size() == 9, "Spin vector should contain one spin per lattice site");

    for (int spin : model.spins)
    {
        passed &= check(spin == 1, "Cold up configuration should contain only +1 spins");
    }

    model.calculate_total_energy();
    model.calculate_total_magnetization();
    passed &= check(nearly_equal(model.Total_energy, -18.0), "3x3 aligned Ising configuration should have total energy -18");
    passed &= check(nearly_equal(model.Total_magnetization, 9.0), "3x3 all-up configuration should have magnetization +9");

    // ---------------------------------------------------------------------
    // Cold all-down configuration
    // ---------------------------------------------------------------------
    model.discrete_spin_generator(generator, false, "down");
    model.calculate_total_energy();
    model.calculate_total_magnetization();
    passed &= check(nearly_equal(model.Total_energy, -18.0), "3x3 all-down configuration should also have total energy -18");
    passed &= check(nearly_equal(model.Total_magnetization, -9.0), "3x3 all-down configuration should have magnetization -9");

    // ---------------------------------------------------------------------
    // Single flipped spin
    //
    // Flipping one spin breaks four aligned bonds.
    // Each bond changes from -1 to +1, giving an energy increase of 8.
    // E = -18 + 8 = -10.
    // ---------------------------------------------------------------------
    model.discrete_spin_generator(generator, false, "up");
    model.spins[0] = -1;
    model.calculate_total_energy();
    model.calculate_total_magnetization();
    passed &= check(nearly_equal(model.Total_energy, -10.0),"Single flipped spin should raise total energy from -18 to -10");
    passed &= check(nearly_equal(model.Total_magnetization, 7.0), "Single flipped spin should reduce magnetization from 9 to 7");

    // ---------------------------------------------------------------------
    // Hot configuration validity
    // ---------------------------------------------------------------------
    model.discrete_spin_generator(generator, true, "hot");
    passed &= check(model.spins.size() == 9, "Hot configuration should contain one spin per lattice site");

    for (int spin : model.spins)
    {
        passed &= check(spin == -1 || spin == 1, "Hot configuration should contain only -1 or +1 spins");
    }

    if (!passed)
    {
        return 1;
    }

    std::cout << "Ising-model unit tests passed.\n";
    return 0;
}
