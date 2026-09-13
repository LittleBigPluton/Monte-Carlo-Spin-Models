# Monte Carlo Spin Models

[![CI](https://github.com/LittleBigPluton/Monte-Carlo-Spin-Models/actions/workflows/ci.yml/badge.svg)](https://github.com/LittleBigPluton/Monte-Carlo-Spin-Models/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)

C++ simulations of two-dimensional lattice spin systems using **Metropolis**, **Wolff cluster updates** and **Hybrid Monte Carlo (HMC)**. The project combines reusable simulation code, reproducible experiment drivers, statistical error analysis, automated unit tests, CI and Python-based visualization.

The repository covers the discrete-spin **2D Ising model** and the continuous-spin **2D XY model**, with emphasis on phase-transition behaviour, Monte Carlo efficiency, autocorrelation and uncertainty estimation.

## Key outputs

| Workflow | Main outputs |
|---|---|
| Ising Metropolis temperature scan | Energy and magnetization observables across temperature |
| Metropolis vs Wolff comparison | Direct comparison of local and cluster updates, including autocorrelation behaviour |
| Error-analysis study | Autocorrelation, blocking and bootstrap uncertainty estimates |
| XY HMC parameter scan | Acceptance rate, autocorrelation times and efficiency proxies across leapfrog settings |
| XY HMC temperature scan | Energy density, magnetization density, standard error, acceptance rate and autocorrelation diagnostics |
| Automated validation | Unit tests for lattice geometry, Ising observables, statistics and core XY-HMC calculations |

### Example HMC outputs

| Parameter scan | Temperature scan |
|---|---|
| ![HMC acceptance-rate parameter scan](plots/hmc_xy/parameter_scan/hmc_xy_parameter_scan_AcceptanceRate.png) | ![XY magnetization density vs temperature](plots/hmc_xy/temperature_scan/hmc_xy_magnetization_density_temperature_scan.png) |
| Acceptance rate across leapfrog settings | Magnetization-density evolution across temperature |

Additional figures are available under [`plots/`](plots/).

## Algorithms

### 2D Ising model

**Metropolis single-spin updates** provide the baseline local Markov-chain Monte Carlo method. Proposed spin flips are accepted according to the Boltzmann weight.

**Wolff single-cluster updates** provide a non-local alternative that can reduce critical slowing down near the Ising critical region by updating correlated clusters rather than individual spins.

### 2D XY model

**Hybrid Monte Carlo** evolves continuous spin angles through leapfrog integration and applies a Metropolis accept/reject step to the resulting trajectory.

The XY workflow includes:

- leapfrog step-size and step-count scans,
- acceptance-rate diagnostics,
- integrated autocorrelation-time estimates,
- temperature scans,
- and efficiency proxies based on autocorrelation and trajectory cost.

## Features

- 2D periodic square lattices
- Cold and hot Ising-spin initialization
- Ising energy and magnetization calculations
- Metropolis single-spin sampling
- Wolff single-cluster sampling
- Classical XY model with Hybrid Monte Carlo
- Leapfrog integration and Hamiltonian evaluation
- Autocorrelation-time estimation
- Blocking and bootstrap error analysis
- Specific-heat and magnetic-susceptibility utilities
- CSV and `.dat` result export
- Python plotting workflows
- Automated unit tests with CTest
- GCC and Clang validation through GitHub Actions

## Repository layout

```text
.
├── .github/
│   └── workflows/
│       └── ci.yml                  # GCC/Clang build and CTest workflow
├── include/                        # Public C++ headers
├── src/                            # Core library implementation
├── experiments/                    # Simulation and analysis drivers
│   ├── metropolis.cpp
│   ├── temperature_scan.cpp
│   ├── wolff_cluster.cpp
│   ├── metropolis_wolff_comparison.cpp
│   ├── error_analysis.cpp
│   ├── hmc_xy_parameter_scan.cpp
│   └── hmc_xy_temperature_scan.cpp
├── tests/
│   └── unit/                       # Automated unit tests
│       ├── test_lattice.cpp
│       ├── test_ising_model.cpp
│       ├── test_statistics.cpp
│       └── test_hmc_xy_model.cpp
├── scripts/                        # Python plotting scripts
├── plots/                          # Generated example figures
├── CMakeLists.txt
├── requirements.txt
├── LICENSE
└── README.md
```

The separation is intentional:

- `src/` and `include/` contain reusable simulation code.
- `experiments/` contains executable scientific workflows and parameter scans.
- `tests/unit/` contains deterministic automated validation used by CTest and CI.

## Core library

The reusable C++ implementation is built as the static library `ising_model_core` and exposed in CMake as `MonteCarlo::Core`.

| Component | Purpose |
|---|---|
| `lattice` | Periodic square-lattice geometry and nearest-neighbour construction |
| `ising_model` | Discrete-spin (`±1`) Ising model |
| `metropolis` | Single-spin Metropolis updates |
| `wolff_cluster` | Wolff single-cluster updates |
| `hmc_xy_model` | Continuous-spin XY model, Hamiltonian evaluation and leapfrog integration |
| `measurements` | Specific heat and magnetic susceptibility from sampled observables |
| `statistics` | Mean, variance and power-law fitting |
| `error_analysis` | Autocorrelation, blocking and bootstrap uncertainty estimation |
| `data_export` | CSV and `.dat` export helpers |

## Requirements

### C++ build

- CMake **3.20 or newer**
- A C++17 compiler
- GCC or Clang recommended on Linux
- OpenMP optional

### Python plotting

- Python 3
- dependencies listed in `requirements.txt`

## Build

Clone the current repository:

```bash
git clone https://github.com/LittleBigPluton/Monte-Carlo-Spin-Models.git
cd Monte-Carlo-Spin-Models
```

Configure and build the library, experiments and unit tests:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DMC_BUILD_TESTS=ON \
    -DMC_BUILD_EXPERIMENTS=ON

cmake --build build --parallel
```

### Useful CMake options

| Option | Default | Purpose |
|---|---:|---|
| `MC_BUILD_TESTS` | `ON` | Build automated unit tests |
| `MC_BUILD_EXPERIMENTS` | `ON` | Build simulation and analysis executables |
| `MC_ENABLE_OPENMP` | `OFF` | Enable OpenMP support |
| `MC_ENABLE_LTO` | `ON` | Enable link-time optimization where supported |
| `MC_ENABLE_NATIVE_OPTIMIZATIONS` | `OFF` | Enable architecture-specific optimization |
| `MC_ENABLE_SANITIZERS` | `OFF` | Enable AddressSanitizer and UndefinedBehaviorSanitizer with GCC/Clang |

For example, a debug build with sanitizers can be configured with:

```bash
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Debug \
    -DMC_BUILD_TESTS=ON \
    -DMC_BUILD_EXPERIMENTS=OFF \
    -DMC_ENABLE_SANITIZERS=ON
```

## Automated tests

The repository contains deterministic unit tests for the reusable simulation components.

```bash
ctest --test-dir build --output-on-failure
```

The current automated test suite covers:

| Test | Validation |
|---|---|
| `test_lattice` | Periodic neighbour construction and lattice indexing |
| `test_ising_model` | Spin initialization, total energy and total magnetization |
| `test_statistics` | Mean, variance, filtering and power-law fitting |
| `test_hmc_xy_model` | XY energy, magnetization, gradient, kinetic energy, Hamiltonian and basic leapfrog behaviour |

The experiment programs are intentionally **not** registered as CTest tests because they are scientific workflows rather than deterministic unit tests.

## Continuous integration

GitHub Actions builds the repository with both **GCC** and **Clang**, compiles the core library and experiment executables and runs the registered CTest suite.

This checks that the project can be configured and built from a clean environment and that the automated validation passes independently of the local development setup.

## Running experiments

Experiment executables are built into `build/`.

### Ising Metropolis simulation

```bash
./build/metropolis
```

### Ising temperature scan

```bash
./build/temperature_scan
```

### Wolff cluster simulation

```bash
./build/wolff_cluster
```

### Metropolis vs Wolff comparison

```bash
./build/metropolis_wolff_comparison
```

### Error-analysis workflow

```bash
./build/error_analysis
```

### XY HMC parameter scan

```bash
./build/hmc_xy_parameter_scan
```

The parameter scan evaluates the built-in leapfrog grid:

| Quantity | Values |
|---|---|
| Leapfrog step size `Δt` | `0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14` |
| Leapfrog steps | `4, 8, 12, 16, 20, 24, 32` |

### XY HMC temperature scan

```bash
./build/hmc_xy_temperature_scan
```

The temperature-scan workflow prompts for:

- lattice size,
- number of equilibration trajectories,
- number of measurement trajectories,
- leapfrog step size,
- and number of leapfrog steps.

The current driver scans temperatures from `0.1` to `3.0`.

## Recommended simulation scales

These are practical starting points for the repository rather than hard-coded requirements.

| Workflow | Lattice size | Equilibration | Measurement length | Notes |
|---|---:|---:|---:|---|
| Ising Metropolis temperature scan | 64 | optional / discarded through filtering | `100000` sweeps | Baseline local-update scan |
| Metropolis vs Wolff comparison | 64 | optional / discarded through filtering | `100000` updates | Compare local and cluster sampling |
| Coarse XY HMC parameter scan | 16 | none explicit in driver | 500–1000 trajectories per point | Identify promising leapfrog settings |
| Final XY HMC parameter scan | 16 | none explicit in driver | 3000–5000 trajectories per point | Refine the parameter choice |
| Coarse XY HMC temperature scan | 16 | 500–1000 trajectories | 2000–3000 trajectories | Use a suitable parameter-scan point |
| Final XY HMC temperature scan | 16 | 1000–3000 trajectories | 5000–10000 trajectories | Higher-statistics temperature study |

## Output files

Simulation outputs are generated at runtime and are not required for building or testing the project.

### XY HMC temperature scan

Summary:

```text
results/hmc_xy/summary/hmc_xy_temperature_scan.csv
```

Per-temperature chains:

```text
results/hmc_xy/chains/energies_<k>.dat
results/hmc_xy/chains/magnetizations_<k>.dat
```

The summary contains:

| Column | Meaning |
|---|---|
| `Temperature` | Simulation temperature |
| `LeapfrogStepSize` | Leapfrog step size `Δt` |
| `LeapfrogSteps` | Leapfrog steps per HMC trajectory |
| `AcceptanceRate` | Fraction of accepted HMC trajectories |
| `EnergyDensity` | Mean energy density |
| `MagnetizationDensity` | Mean magnetization density |
| `MagnetizationSE` | Autocorrelation-corrected standard error |
| `EnergyTau` | Integrated autocorrelation time for energy density |
| `MagnetizationTau` | Integrated autocorrelation time for magnetization density |

### XY HMC parameter scan

Summary:

```text
results/hmc_xy/summary/hmc_xy_parameter_scan.csv
```

The parameter-scan output includes:

| Column | Meaning |
|---|---|
| `Temperature` | Fixed scan temperature |
| `LeapfrogStepSize` | Leapfrog step size `Δt` |
| `LeapfrogSteps` | Number of leapfrog steps |
| `AcceptanceRate` | Fraction of accepted trajectories |
| `EnergyTau` | Integrated autocorrelation time for energy density |
| `MagnetizationTau` | Integrated autocorrelation time for magnetization density |
| `IndependentTimeEnergy` | Cost proxy using energy autocorrelation and trajectory length |
| `IndependentTimeMagnetization` | Cost proxy using magnetization autocorrelation and trajectory length |
| `MagnetizationDensity` | Mean magnetization density |
| `MagnetizationSE` | Autocorrelation-corrected standard error |

Other Ising workflows generate outputs under directories such as:

```text
results/summary/
results/chains/
results/autocorrelation/
results/comparison/
```

## Interpreting HMC diagnostics

- Higher `AcceptanceRate` generally indicates more stable numerical integration, although acceptance alone does not determine sampling efficiency.
- Lower `EnergyTau` and `MagnetizationTau` indicate faster decorrelation of the corresponding observables.
- Lower `IndependentTimeEnergy` and `IndependentTimeMagnetization` indicate a lower estimated computational cost per effectively independent sample.
- `MagnetizationDensity` together with `MagnetizationSE` provides the main order/disorder observable reported by the XY temperature scan.

## Python plotting

Create an isolated environment:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

After generating simulation data, the HMC figures can be reproduced with:

```bash
python3 scripts/plot_hmc_xy_parameter_scan.py
python3 scripts/plot_hmc_xy_temperature_scan.py
```

Generated HMC figures are stored under:

```text
plots/hmc_xy/parameter_scan/
plots/hmc_xy/temperature_scan/
```

Examples include:

- acceptance rate vs leapfrog step size,
- energy and magnetization autocorrelation times,
- independent-time efficiency diagnostics,
- magnetization density vs temperature,
- energy density vs temperature,
- and combined temperature-scan diagnostics.

## Scope

This repository is intended as a computational-statistical-physics project and numerical Monte Carlo implementation. It emphasizes transparent algorithms, reproducible experiment drivers, statistical diagnostics and automated validation.

The included workflows are designed for numerical experimentation and comparison of Monte Carlo methods rather than as a general-purpose high-performance simulation framework.

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for details.
