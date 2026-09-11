#include<iostream>
#include<vector>
#include<random>
#include<chrono>
#include<fstream>
#include<filesystem>

#include"hmc_xy_model.hpp"
#include"statistics.hpp"
#include"error_analysis.hpp"

class HMCXYParameterScan{
private:
  HMCXYModel model;
  std::mt19937 generator;
  int trajectories = 0;
  double filter = 0.1;

public:
  void scan(){
    std::filesystem::create_directories("results/hmc_xy/summary");

    std::ofstream file("results/hmc_xy/summary/hmc_xy_parameter_scan.csv",std::ios::out);
    file<<"Temperature,LeapfrogStepSize,LeapfrogSteps,AcceptanceRate,EnergyTau,MagnetizationTau,IndependentTimeEnergy,IndependentTimeMagnetization,MagnetizationDensity,MagnetizationSE"<<std::endl;

    model.initialize();

    unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
    generator.seed(seed);

    std::cout<<"Enter the number of HMC trajectories for parameter scan: "<<std::endl;
    std::cin>>trajectories;

    if(trajectories<=0){std::cerr<<"Number of trajectories must be positive."<<std::endl; return;}

    double temperature = 1.0;
    std::cout<<"Enter temperature for parameter scan: "<<std::endl;
    std::cin>>temperature;

    std::vector<double> leapfrog_step_sizes = {0.02,0.04,0.06,0.08,0.10,0.12,0.14};
    std::vector<int> leapfrog_step_counts = {4,8,12,16,20,24,32};

    for(double leapfrog_step_size:leapfrog_step_sizes){
      for(int leapfrog_steps:leapfrog_step_counts){
        std::cout<<"Temperature = "<<temperature<<", step size = "<<leapfrog_step_size<<", leapfrog steps = "<<leapfrog_steps<<std::endl;

        model.xy_spins.clear();
        model.initialize_xy_spins(generator);

        std::vector<double> energies;
        std::vector<double> magnetizations;
        energies.reserve(trajectories);
        magnetizations.reserve(trajectories);

        int accepted_trajectories = 0;

        for(int trajectory=0;trajectory<trajectories;trajectory++){
          bool accepted = model.hmc_trajectory(temperature,leapfrog_steps,leapfrog_step_size,generator);

          if(accepted){accepted_trajectories++;}

          energies.push_back(model.energy_density());
          magnetizations.push_back(model.magnetization_density());
        }

        double acceptance_rate = static_cast<double>(accepted_trajectories)/trajectories;

        std::vector<double> decay_values;
        double energy_tau = autocorrelation_time(energies,decay_values,filter);
        decay_values.clear();

        double magnetization_tau = autocorrelation_time(magnetizations,decay_values,filter);
        decay_values.clear();

        double magnetization_average = mean_value(magnetizations,filter);
        double magnetization_variance = variance(magnetizations,magnetization_average,filter);
        double magnetization_standard_error = ac_standard_error(magnetization_tau,magnetization_variance,static_cast<int>(magnetizations.size()),filter);

        double independent_time_energy = energy_tau*leapfrog_steps;
        double independent_time_magnetization = magnetization_tau*leapfrog_steps;

        file<<temperature<<","<<leapfrog_step_size<<","<<leapfrog_steps<<","<<acceptance_rate<<","
            <<energy_tau<<","<<magnetization_tau<<","<<independent_time_energy<<","
            <<independent_time_magnetization<<","<<magnetization_average<<","<<magnetization_standard_error<<std::endl;
      }
    }

    file.close();
  }
};

int main(){
  HMCXYParameterScan scan;
  scan.scan();
  return 0;
}
