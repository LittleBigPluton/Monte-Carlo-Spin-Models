#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<chrono>
#include<fstream>
#include<filesystem>
#include<cmath>

#include"ising_model.hpp"
#include"wolff_cluster.hpp"
#include"measurements.hpp"
#include"statistics.hpp"
#include"error_analysis.hpp"
#include"data_export.hpp"

class WolffPhaseScan{
private:
  Ising model;
  WolffCluster wolff;
  Measurements measurements;
  MeasurementSummary summary;
  std::vector<double> cluster_size_densities;
  std::mt19937 generator;
  int cluster_flips = 0;
  double filter = 0.1;

  void error_calculations(std::vector<double>& measurement,double measurement_variance,double beta,double& tau,double& standard_error,double& secondary_error){
    std::vector<double> correlation_function;
    tau = autocorrelation_time(measurement,correlation_function,filter);
    standard_error = ac_standard_error(tau,measurement_variance,cluster_flips,filter)/model.grid_points;
    std::vector<double> effective_correlation_function;
    secondary_error = error_propagation(measurement,effective_correlation_function,filter)/model.grid_points/beta;
  }

public:
  void phase_region(){
    std::filesystem::create_directories("results/wolff/summary");
    std::filesystem::create_directories("results/wolff/chains");
    std::filesystem::create_directories("results/wolff/cluster");

    /****************************************************************/
    /*        Evaluate phase region with Wolff cluster algorithm     */
    /****************************************************************/

    //Create a file to save observables data for phase scan
    std::ofstream phase_scan_cluster("results/wolff/summary/phase_scan_cluster.dat",std::ios::out);
    phase_scan_cluster<<"Temperature AED SH AMD MS"<<std::endl;

    //Create a file to save errors and taus for phase scan
    std::ofstream stats_phase_scan_cl("results/wolff/summary/stats_phase_scan_cl.dat",std::ios::out);
    stats_phase_scan_cl<<"Temperature E_Tau E_SE SH_E M_Tau M_SE MS_E CLDav CLDvar"<<std::endl;

    //Get lattice size value from user
    model.initialize();

    //Create a unique seed for random numbers
    unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
    generator.seed(seed);

    //Get number of cluster flips from user
    std::cout<<"Enter the number of cluster flips: "<<std::endl;
    std::cin>>cluster_flips;

    if(cluster_flips<=0){std::cerr<<"Number of cluster flips must be positive."<<std::endl; return;}

    for(int temp_index=210;temp_index<=240;temp_index+=2){
      double Temp = temp_index/100.0;

      //Print current temperature
      std::cout<<"Temperature = "<<Temp<<std::endl;

      //Clear vectors before each temperature
      measurements.energies.clear();
      measurements.magnetizations.clear();
      cluster_size_densities.clear();
      summary = MeasurementSummary{};

      //Create initial spin configuration
      model.discrete_spin_generator(generator,false,"cold");
      model.calculate_total_energy();
      model.calculate_total_magnetization();

      for(int flip=0;flip<cluster_flips;flip++){
        //Perform one Wolff cluster flip
        wolff.sweep(model,cluster_size_densities,Temp,generator);

        //Recalculate total energy after one cluster flip and store it
        model.calculate_total_energy();
        measurements.energies.push_back(model.Total_energy);

        //Recalculate total magnetization after one cluster flip and store it
        model.calculate_total_magnetization();
        measurements.magnetizations.push_back(std::abs(model.Total_magnetization));
      }
      // Export chains and cluster density
      std::string temp_tag = std::to_string(static_cast<int>(Temp*100));
      dataexport("results/wolff/chains/energies_"+temp_tag,"dat",measurements.energies,"ClusterFlip","Energies",filter,100,model.grid_points,true);
      dataexport("results/wolff/chains/magnetizations_"+temp_tag,"dat",measurements.magnetizations,"ClusterFlip","Magnetizations",filter,100,model.grid_points,true);
      dataexport("results/wolff/cluster/cluster_density_"+temp_tag,"dat",cluster_size_densities,"ClusterFlip","ClusterDensity",0.0,1,model.grid_points,false);

      //Calculate primary observables
      double en_av = mean_value(measurements.energies,filter);
      double en_var = variance(measurements.energies,en_av,filter);
      double mag_av = mean_value(measurements.magnetizations,filter);
      double mag_var = variance(measurements.magnetizations,mag_av,filter);

      summary.energy_density = en_av/model.grid_points;
      summary.magnetization_density = mag_av/model.grid_points;

      //Calculate secondary observables
      calculate_specific_heat(en_var,Temp,model.grid_points,summary);
      calculate_magnetic_susceptibility(mag_var,Temp,model.grid_points,summary);

      phase_scan_cluster<<Temp<<" "<<summary.energy_density<<" "<<summary.specific_heat<<" "
                        <<summary.magnetization_density<<" "<<summary.magnetic_susceptibility<<std::endl;

      //Calculate cluster-size density average and variance
      double cl_den_av = mean_value(cluster_size_densities,0.0);
      double cl_den_var = variance(cluster_size_densities,cl_den_av,0.0);

      //Calculate integrated autocorrelation time and standard errors
      double E_tau = 0.0, M_tau = 0.0, E_SE = 0.0, M_SE = 0.0, SH_E = 0.0, MS_E = 0.0;
      error_calculations(measurements.energies,en_var,Temp*Temp,E_tau,E_SE,SH_E);
      error_calculations(measurements.magnetizations,mag_var,Temp,M_tau,M_SE,MS_E);

      stats_phase_scan_cl<<Temp<<" "<<E_tau<<" "<<E_SE<<" "<<SH_E<<" "
                         <<M_tau<<" "<<M_SE<<" "<<MS_E<<" "<<cl_den_av<<" "<<cl_den_var<<std::endl;
    }

    phase_scan_cluster.close();
    stats_phase_scan_cl.close();
  }
};

int main(){
  WolffPhaseScan obj;
  obj.phase_region();
  return 0;
}
