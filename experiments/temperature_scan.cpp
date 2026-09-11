#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<chrono>
#include<fstream>
#include<filesystem>
#include<cmath>

#include"ising_model.hpp"
#include"metropolis.hpp"
#include"measurements.hpp"
#include"statistics.hpp"
#include"error_analysis.hpp"
#include"data_export.hpp"

class TemperatureScan{
private:
  Ising model;
  Metropolis metropolis;
  Measurements measurements;
  MeasurementSummary summary;
  std::mt19937 generator;
  int Nmc = 0;
  double filter = 0.1;

  void create_directories(){
    std::filesystem::create_directories("results/temperature_scan/summary");
    std::filesystem::create_directories("results/temperature_scan/chains");
  }

  void clear_measurements(){
    measurements.energies.clear();
    measurements.magnetizations.clear();
    summary = MeasurementSummary{};
  }

  void recalculate_totals(){
    model.Total_energy = 0;
    model.Total_magnetization = 0;
    model.calculate_total_energy();
    model.calculate_total_magnetization();
  }

  void clean_start(){
    clear_measurements();
    model.spins.clear();
    model.discrete_spin_generator(generator,false,"cold");
    recalculate_totals();
  }

  void calculate_observables(double& Temp,double& energy_variance,double& magnetization_variance){
    double energy_mean = mean_value(measurements.energies,filter);
    energy_variance = variance(measurements.energies,energy_mean,filter);
    double magnetization_mean = mean_value(measurements.magnetizations,filter);
    magnetization_variance = variance(measurements.magnetizations,magnetization_mean,filter);
    summary.energy_density = energy_mean/model.grid_points;
    summary.magnetization_density = magnetization_mean/model.grid_points;
    calculate_specific_heat(energy_variance,Temp,model.grid_points,summary);
    calculate_magnetic_susceptibility(magnetization_variance,Temp,model.grid_points,summary);
  }

  void error_calculations(std::vector<double>& measurement,double measurement_variance,double secondary_scale,double& tau,double& standard_error,double& secondary_error){
    std::vector<double> correlation_function;
    tau = autocorrelation_time(measurement,correlation_function,filter);
    standard_error = ac_standard_error(tau,measurement_variance,Nmc,filter)/model.grid_points;
    std::vector<double> effective_correlation_function;
    secondary_error = error_propagation(measurement,effective_correlation_function,filter)/model.grid_points/secondary_scale;
  }

public:
  void temperature_scan(){
    create_directories();

    //Create a file to save observables data
    std::ofstream observables_scan("results/temperature_scan/summary/observables_scan.dat",std::ios::out);
    observables_scan<<"Temperature AED SH AMD MS"<<std::endl;

    //Create a file to save errors and taus for temperature scan
    std::ofstream stats_temp_scan("results/temperature_scan/summary/stats_temp_scan.dat",std::ios::out);
    stats_temp_scan<<"Temperature E_Tau E_SE SH_E M_Tau M_SE MS_E"<<std::endl;

    //Get lattice size value from user
    model.initialize();

    //Create a unique seed for random numbers
    unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
    generator.seed(seed);

    //Get Monte Carlo chain length
    std::cout<<"Enter the number of Markov sweeps for temperature scan = "<<std::endl;
    std::cin>>Nmc;

    if(Nmc<=0){std::cerr<<"Number of Markov sweeps must be positive."<<std::endl; return;}

    for(int temp_index=10;temp_index<=40;temp_index++){
      double Temp = temp_index/10.0;

      //Print current temperature
      std::cout<<"Temperature = "<<Temp<<std::endl;

      clean_start();
      metropolis.lookup(Temp);

      auto start = std::chrono::high_resolution_clock::now();

      for(int sweep=0;sweep<Nmc;sweep++){
        metropolis.sweep(model,generator,measurements);
      }

      //Calculate statistical values for energies and magnetizations
      double energy_variance = 0.0;
      double magnetization_variance = 0.0;
      calculate_observables(Temp,energy_variance,magnetization_variance);

      observables_scan<<Temp<<" "<<summary.energy_density<<" "<<summary.specific_heat<<" "<<summary.magnetization_density<<" "<<summary.magnetic_susceptibility<<std::endl;

      //Calculate autocorrelation time, standard error, and propagated secondary-observable errors
      double E_tau = 0.0, E_SE = 0.0, SH_E = 0.0;
      double M_tau = 0.0, M_SE = 0.0, MS_E = 0.0;

      error_calculations(measurements.energies,energy_variance,Temp*Temp,E_tau,E_SE,SH_E);
      error_calculations(measurements.magnetizations,magnetization_variance,Temp,M_tau,M_SE,MS_E);

      stats_temp_scan<<Temp<<" "<<E_tau<<" "<<E_SE<<" "<<SH_E<<" "<<M_tau<<" "<<M_SE<<" "<<MS_E<<std::endl;

      std::string temp_tag = std::to_string(temp_index);
      dataexport("results/temperature_scan/chains/energies_"+temp_tag,"dat",measurements.energies,"MonteCarloTime","Energies",filter,100,model.grid_points,true);
      dataexport("results/temperature_scan/chains/magnetizations_"+temp_tag,"dat",measurements.magnetizations,"MonteCarloTime","Magnetizations",filter,100,model.grid_points,true);

      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
      std::cout<<"Temperature scan time: "<<duration.count()<<" microseconds"<<std::endl;
    }

    observables_scan.close();
    stats_temp_scan.close();
  }
};

int main(){
  TemperatureScan scan;
  scan.temperature_scan();
  return 0;
}
