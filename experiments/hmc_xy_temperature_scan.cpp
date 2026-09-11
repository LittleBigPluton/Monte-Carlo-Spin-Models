#include<iostream>
#include<vector>
#include<random>
#include<chrono>
#include<fstream>
#include<filesystem>
#include<string>

#include"hmc_xy_model.hpp"
#include"statistics.hpp"
#include"error_analysis.hpp"
#include"data_export.hpp"

class HMCXYTemperatureScan{
private:
  HMCXYModel model;
  std::mt19937 generator;
  int trajectories = 0;
  int equilibration_steps = 0;
  int leapfrog_steps = 0;
  double leapfrog_step_size = 0.0;
  double filter = 0.0;

public:
  void scan(){
    std::filesystem::create_directories("results/hmc_xy/summary");
    std::filesystem::create_directories("results/hmc_xy/chains");

    std::ofstream file("results/hmc_xy/summary/hmc_xy_temperature_scan.csv",std::ios::out);
    file<<"Temperature,LeapfrogStepSize,LeapfrogSteps,AcceptanceRate,EnergyDensity,MagnetizationDensity,MagnetizationSE,EnergyTau,MagnetizationTau"<<std::endl;

    model.initialize();

    unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
    generator.seed(seed);

    std::cout<<"Enter the number of equilibration trajectories: "<<std::endl;
    std::cin>>equilibration_steps;

    std::cout<<"Enter the number of measurement trajectories: "<<std::endl;
    std::cin>>trajectories;

    std::cout<<"Enter leapfrog step size: "<<std::endl;
    std::cin>>leapfrog_step_size;

    std::cout<<"Enter number of leapfrog steps: "<<std::endl;
    std::cin>>leapfrog_steps;

    if(equilibration_steps<0){std::cerr<<"Equilibration steps cannot be negative."<<std::endl; return;}
    if(trajectories<=0){std::cerr<<"Number of trajectories must be positive."<<std::endl; return;}
    if(leapfrog_step_size<=0.0){std::cerr<<"Leapfrog step size must be positive."<<std::endl; return;}
    if(leapfrog_steps<=0){std::cerr<<"Number of leapfrog steps must be positive."<<std::endl; return;}

    for(int temp_index=1;temp_index<=30;temp_index++){
      double temperature = temp_index/10.0;

      std::cout<<"Temperature = "<<temperature<<std::endl;

      model.xy_spins.clear();
      model.initialize_xy_spins(generator);

      for(int step=0;step<equilibration_steps;step++){
        model.hmc_trajectory(temperature,leapfrog_steps,leapfrog_step_size,generator);
      }

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

      double energy_average = mean_value(energies,filter);
      double magnetization_average = mean_value(magnetizations,filter);
      double magnetization_variance = variance(magnetizations,magnetization_average,filter);

      std::vector<double> decay_values;
      double energy_tau = autocorrelation_time(energies,decay_values,filter);
      decay_values.clear();

      double magnetization_tau = autocorrelation_time(magnetizations,decay_values,filter);
      decay_values.clear();

      double magnetization_standard_error = ac_standard_error(magnetization_tau,magnetization_variance,static_cast<int>(magnetizations.size()),filter);

      file<<temperature<<","<<leapfrog_step_size<<","<<leapfrog_steps<<","<<acceptance_rate<<","
          <<energy_average<<","<<magnetization_average<<","<<magnetization_standard_error<<","
          <<energy_tau<<","<<magnetization_tau<<std::endl;

      std::string temp_tag = std::to_string(temp_index);
      dataexport("results/hmc_xy/chains/energies_"+temp_tag,"dat",energies,"Trajectory","EnergyDensity",0.0,100,model.grid_points,false);
      dataexport("results/hmc_xy/chains/magnetizations_"+temp_tag,"dat",magnetizations,"Trajectory","MagnetizationDensity",0.0,100,model.grid_points,false);
    }

    file.close();
  }
};

int main(){
  HMCXYTemperatureScan scan;
  scan.scan();
  return 0;
}
