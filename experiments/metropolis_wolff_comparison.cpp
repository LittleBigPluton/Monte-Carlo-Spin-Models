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
#include"wolff_cluster.hpp"
#include"measurements.hpp"
#include"statistics.hpp"
#include"error_analysis.hpp"
#include"data_export.hpp"

class MetropolisWolffComparison{
private:
  Ising model;
  Metropolis metropolis;
  WolffCluster wolff;
  Measurements measurements;
  MeasurementSummary summary;
  std::vector<double> cluster_size_densities;
  std::mt19937 generator;
  int Nmc = 0;
  int cluster_flips = 0;
  double filter = 0.1;

  void create_directories(){
    std::filesystem::create_directories("results/comparison/summary");
    std::filesystem::create_directories("results/comparison/metropolis/chains");
    std::filesystem::create_directories("results/comparison/wolff/chains");
    std::filesystem::create_directories("results/comparison/wolff/cluster");
  }

  void clear_measurements(){
    measurements.energies.clear();
    measurements.magnetizations.clear();
    cluster_size_densities.clear();
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

  void error_calculations(std::vector<double>& measurement,double measurement_variance,int sample_count,double secondary_scale,double& tau,double& standard_error,double& secondary_error){
    std::vector<double> correlation_function;
    tau = autocorrelation_time(measurement,correlation_function,filter);
    standard_error = ac_standard_error(tau,measurement_variance,sample_count,filter)/model.grid_points;
    std::vector<double> effective_correlation_function;
    secondary_error = error_propagation(measurement,effective_correlation_function,filter)/model.grid_points/secondary_scale;
  }

public:
  void run(){
    create_directories();
    model.initialize();

    unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
    generator.seed(seed);

    std::cout<<"Enter the number of Wolff cluster flips: "<<std::endl;
    std::cin>>cluster_flips;
    if(cluster_flips<=0){std::cerr<<"Number of Wolff cluster flips must be positive."<<std::endl; return;}

    std::cout<<"Enter the number of Metropolis sweeps: "<<std::endl;
    std::cin>>Nmc;
    if(Nmc<=0){std::cerr<<"Number of Metropolis sweeps must be positive."<<std::endl; return;}

    run_wolff_phase_scan();
    run_metropolis_phase_scan();
  }

  void run_wolff_phase_scan(){
    std::ofstream phase_scan_cluster("results/comparison/summary/phase_scan_cluster.dat",std::ios::out);
    phase_scan_cluster<<"Temperature AED SH AMD MS"<<std::endl;

    std::ofstream stats_phase_scan_cl("results/comparison/summary/stats_phase_scan_cl.dat",std::ios::out);
    stats_phase_scan_cl<<"Temperature E_Tau E_SE SH_E M_Tau M_SE MS_E CLDav CLDvar"<<std::endl;

    for(int temp_index=210;temp_index<=240;temp_index+=2){
      double Temp = temp_index/100.0;
      std::cout<<"Wolff temperature = "<<Temp<<std::endl;

      clean_start();

      auto start = std::chrono::high_resolution_clock::now();

      for(int flip=0;flip<cluster_flips;flip++){
        wolff.sweep(model,cluster_size_densities,Temp,generator);
        recalculate_totals();
        measurements.energies.push_back(model.Total_energy);
        measurements.magnetizations.push_back(std::abs(model.Total_magnetization));
      }

      double energy_variance = 0.0;
      double magnetization_variance = 0.0;
      calculate_observables(Temp,energy_variance,magnetization_variance);

      phase_scan_cluster<<Temp<<" "<<summary.energy_density<<" "<<summary.specific_heat<<" "<<summary.magnetization_density<<" "<<summary.magnetic_susceptibility<<std::endl;

      double cluster_density_average = mean_value(cluster_size_densities,0.0);
      double cluster_density_variance = variance(cluster_size_densities,cluster_density_average,0.0);

      double E_tau = 0.0, E_SE = 0.0, SH_E = 0.0;
      double M_tau = 0.0, M_SE = 0.0, MS_E = 0.0;

      error_calculations(measurements.energies,energy_variance,cluster_flips,Temp*Temp,E_tau,E_SE,SH_E);
      error_calculations(measurements.magnetizations,magnetization_variance,cluster_flips,Temp,M_tau,M_SE,MS_E);

      stats_phase_scan_cl<<Temp<<" "<<E_tau<<" "<<E_SE<<" "<<SH_E<<" "<<M_tau<<" "<<M_SE<<" "<<MS_E<<" "<<cluster_density_average<<" "<<cluster_density_variance<<std::endl;

      std::string temp_tag = std::to_string(temp_index);
      dataexport("results/comparison/wolff/chains/energies_"+temp_tag,"dat",measurements.energies,"ClusterFlip","Energies",filter,100,model.grid_points,true);
      dataexport("results/comparison/wolff/chains/magnetizations_"+temp_tag,"dat",measurements.magnetizations,"ClusterFlip","Magnetizations",filter,100,model.grid_points,true);
      dataexport("results/comparison/wolff/cluster/cluster_density_"+temp_tag,"dat",cluster_size_densities,"ClusterFlip","ClusterDensity",0.0,100,model.grid_points,false);

      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
      std::cout<<"Wolff time: "<<duration.count()<<" microseconds"<<std::endl;
    }

    phase_scan_cluster.close();
    stats_phase_scan_cl.close();
  }

  void run_metropolis_phase_scan(){
    std::ofstream phase_scan_mtrpls("results/comparison/summary/phase_scan_mtrpls.dat",std::ios::out);
    phase_scan_mtrpls<<"Temperature AED SH AMD MS"<<std::endl;

    std::ofstream stats_phase_scan_mp("results/comparison/summary/stats_phase_scan_mp.dat",std::ios::out);
    stats_phase_scan_mp<<"Temperature E_Tau E_SE SH_E M_Tau M_SE MS_E"<<std::endl;

    for(int temp_index=210;temp_index<=240;temp_index+=2){
      double Temp = temp_index/100.0;
      std::cout<<"Metropolis temperature = "<<Temp<<std::endl;

      clean_start();
      metropolis.lookup(Temp);

      auto start = std::chrono::high_resolution_clock::now();

      for(int sweep=0;sweep<Nmc;sweep++){
        metropolis.sweep(model,generator,measurements);
      }

      double energy_variance = 0.0;
      double magnetization_variance = 0.0;
      calculate_observables(Temp,energy_variance,magnetization_variance);

      phase_scan_mtrpls<<Temp<<" "<<summary.energy_density<<" "<<summary.specific_heat<<" "<<summary.magnetization_density<<" "<<summary.magnetic_susceptibility<<std::endl;

      double E_tau = 0.0, E_SE = 0.0, SH_E = 0.0;
      double M_tau = 0.0, M_SE = 0.0, MS_E = 0.0;

      error_calculations(measurements.energies,energy_variance,Nmc,Temp*Temp,E_tau,E_SE,SH_E);
      error_calculations(measurements.magnetizations,magnetization_variance,Nmc,Temp,M_tau,M_SE,MS_E);

      stats_phase_scan_mp<<Temp<<" "<<E_tau<<" "<<E_SE<<" "<<SH_E<<" "<<M_tau<<" "<<M_SE<<" "<<MS_E<<std::endl;

      std::string temp_tag = std::to_string(temp_index);
      dataexport("results/comparison/metropolis/chains/energies_"+temp_tag,"dat",measurements.energies,"MonteCarloTime","Energies",filter,100,model.grid_points,true);
      dataexport("results/comparison/metropolis/chains/magnetizations_"+temp_tag,"dat",measurements.magnetizations,"MonteCarloTime","Magnetizations",filter,100,model.grid_points,true);

      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
      std::cout<<"Metropolis time: "<<duration.count()<<" microseconds"<<std::endl;
    }

    phase_scan_mtrpls.close();
    stats_phase_scan_mp.close();
  }
};

int main(){
  MetropolisWolffComparison comparison;
  comparison.run();
  return 0;
}
