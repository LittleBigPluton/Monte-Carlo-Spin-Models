#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <filesystem>

#include "ising_model.hpp"
#include "metropolis.hpp"
#include "measurements.hpp"
#include "statistics.hpp"
#include "error_analysis.hpp"
#include "data_export.hpp"

class MCS
{
private:
  Ising model;
  Metropolis metropolis;
  Measurements measurements;
  MeasurementSummary summary;
  std::mt19937 generator;

  //Initialize number of bootstrap samples and number of blocks
  int M = 1000, NB = 20;
  //Initialize unequilibrated part length
  double filter = 0.1;
  //Length of Monte Carlo chain
  int Nmc = 0;

public:
  void simulation()
  {
    std::filesystem::create_directories("results/summary");
    std::filesystem::create_directories("results/chains");
    std::filesystem::create_directories("results/autocorrelation");

    //Create a file to hold primary and secondary data values
    std::ofstream observables("results/summary/observables_error_analysis.dat", std::ios::out);
    observables << "Temperature AED AMD SH MS E_tau E_SE M_tau M_SE" << std::endl;

    //Create a file to store errors from error propagation, blocking method, and bootstrap method
    std::ofstream errors("results/summary/errors.dat", std::ios::out);
    errors << "Temperature type EP BM BSM" << std::endl;

    //Initialize lattice and grid points
    model.initialize();

    std::cout << "Enter length of Monte Carlo Chains: ";
    std::cin >> Nmc;

    if(Nmc <= 0){std::cerr << "Number of Monte Carlo sweeps must be positive." << std::endl; return;}

    //Define temperature values
    double Temperature[3] = {2.0, 2.3, 2.6};

    for(double Temp : Temperature)
    {
      //Create lookup table for the given temperature
      metropolis.lookup(Temp);

      //Make separable line on the terminal
      std::cout << std::string(50, '*') << std::endl;

      //Clear everything before use to prevent miscalculations
      measurements.energies.clear();
      measurements.magnetizations.clear();

      //Generate seed
      unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
      generator.seed(seed);

      //Generate spins and calculate initial values
      model.discrete_spin_generator(generator, false, "cold");
      model.calculate_total_energy();
      model.calculate_total_magnetization();

      std::cout << "Temperature = " << Temp << std::endl;

      //Run Metropolis algorithm
      for(int i = 0; i < Nmc; ++i){metropolis.sweep(model, generator, measurements);}

      //Calculate average and variance of the equilibrated part
      double en_av = mean_value(measurements.energies, filter);
      double en_var = variance(measurements.energies, en_av, filter);
      double mag_av = mean_value(measurements.magnetizations, filter);
      double mag_var = variance(measurements.magnetizations, mag_av, filter);

      //Calculate average energy density
      summary.energy_density = en_av / model.grid_points;
      std::cout << "Average Energy Density = " << summary.energy_density << std::endl;

      //Calculate average magnetization density
      summary.magnetization_density = mag_av / model.grid_points;
      std::cout << "Average Magnetization Density = " << summary.magnetization_density << std::endl;

      //Calculate specific heat and magnetic susceptibility
      calculate_specific_heat(en_var, Temp, model.grid_points, summary);
      calculate_magnetic_susceptibility(mag_var, Temp, model.grid_points, summary);

      std::cout << "Specific Heat = " << summary.specific_heat << std::endl
                << "Magnetic Susceptibility = " << summary.magnetic_susceptibility << std::endl;

      observables << Temp << " " << summary.energy_density << " " << summary.magnetization_density << " "
                  << summary.specific_heat << " " << summary.magnetic_susceptibility;

      //Export average energy and magnetization chains
      dataexport("results/chains/energies_" + std::to_string(static_cast<int>(Temp * 10)), "dat",
                 measurements.energies, "Times", "Energies", filter, 100, model.grid_points, true);

      dataexport("results/chains/magnetizations_" + std::to_string(static_cast<int>(Temp * 10)), "dat",
                 measurements.magnetizations, "Times", "Magnetizations", filter, 100, model.grid_points, true);

      //Error calculations for primary and secondary observables
      double E_tau = 0.0, M_tau = 0.0, E_SE = 0.0, M_SE = 0.0;

      error_calculations(measurements.energies, en_var, 'E', "Energy", "SH", Temp * Temp, Temp, E_tau, E_SE);
      observables << " " << E_tau << " " << E_SE;

      error_calculations(measurements.magnetizations, mag_var, 'M', "Magnetization", "MS", Temp, Temp, M_tau, M_SE);
      observables << " " << M_tau << " " << M_SE << std::endl;
    }
  }

  void error_calculations(std::vector<double>& measurement, double measurement_variance, char type,
                          std::string primary, std::string secondary, double Beta, double Temp,
                          double& tau, double& standard_error)
  {
    //Create a vector to store decay values of primary measurement
    std::vector<double> correlation_function;

    //Calculate integrated autocorrelation time
    tau = autocorrelation_time(measurement, correlation_function, filter);

    //Export decay values into a file
    dataexport("results/autocorrelation/" + primary + "_decay_" + std::to_string(static_cast<int>(Temp * 10)), "dat",
               correlation_function, "Times", "Decay", 0.0, 1, model.grid_points, false);

    std::cout << primary << " Tau = " << tau << std::endl;

    //Calculate standard error on the primary measurement; divide by grid points to calculate density error
    standard_error = ac_standard_error(tau, measurement_variance, Nmc, filter) / model.grid_points;
    std::cout << primary << " Standard Error = " << standard_error << std::endl;

    //Create a vector to store effective observable decay values
    std::vector<double> effective_correlation_function;

    //Apply error propagation for secondary observable; divide by grid points and Beta
    double error_prop = error_propagation(measurement, effective_correlation_function, filter) / model.grid_points / Beta;
    std::cout << secondary << " Error propagation = " << error_prop << std::endl;

    //Export effective observable decay values
    dataexport("results/autocorrelation/" + secondary + "_decay_" + std::to_string(static_cast<int>(Temp * 10)), "dat",
               effective_correlation_function, "Time", "Decay", 0.0, 1, model.grid_points, false);

    //Blocking method error on secondary observable
    double blocking_method_variance = blocking_method(measurement, NB, filter) / model.grid_points / Beta;
    std::cout << secondary << " Blocking method = " << blocking_method_variance << std::endl;

    //Bootstrap method error on secondary observable
    double bootstrap_method_variance = bootstrap_method(measurement, M, tau, filter) / model.grid_points / Beta;
    std::cout << secondary << " Bootstrap method = " << bootstrap_method_variance << std::endl;

    //Export all error values into a file
    std::ofstream errors("results/summary/errors.dat", std::ios::app);
    errors << Temp << " " << type << " " << error_prop << " " << blocking_method_variance << " "
           << bootstrap_method_variance << std::endl;
  }
};

int main()
{
  MCS obj;
  obj.simulation();
  return 0;
}
