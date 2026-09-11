#include "lattice.hpp"

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

class TestLattice : public Lattice {
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
    TestLattice lattice;
    lattice.configure(3);
    passed &= check(lattice.grid_points == 9, "3x3 lattice should contain 9 grid points");
    passed &= check(lattice.neighbours.size() == 9, "Neighbour table should contain one entry per lattice site");

    for (std::size_t site = 0; site < lattice.neighbours.size(); ++site)
    {
        passed &= check(lattice.neighbours[site].size() == 4, "Every lattice site should have four nearest neighbours");
        for (int neighbour : lattice.neighbours[site])
        {
            passed &= check(neighbour >= 0 && neighbour < lattice.grid_points, "Neighbour index should remain inside lattice bounds");
        }
    }

    // Site 0 = top-left corner.
    passed &= check(lattice.neighbours[0] == std::vector<int>({2, 6, 1, 3}), "Periodic neighbours of site 0 are incorrect");

    // Site 4 = centre.
    passed &= check(lattice.neighbours[4] == std::vector<int>({3, 1, 5, 7}), "Neighbours of centre site 4 are incorrect");

    // Site 8 = bottom-right corner.
    passed &= check(lattice.neighbours[8] == std::vector<int>({7, 5, 6, 2}), "Periodic neighbours of site 8 are incorrect");

    if (!passed)
    {
        return 1;
    }

    std::cout << "Lattice unit tests passed.\n";
    return 0;
}
