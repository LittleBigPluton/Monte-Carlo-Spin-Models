#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<chrono>
#include<fstream>
#include<filesystem>

#include"ising_model.hpp"
#include"metropolis.hpp"
#include"measurements.hpp"
#include"statistics.hpp"
#include"data_export.hpp"

class MCS{
private:
  Ising model;
  Metropolis metropolis;
  Measurements measurements;
  MeasurementSummary summary;
  std::mt19937 generator;
  int Nmc = 0;
  double filter = 0.1;

public:
  void simulation(){
    std::filesystem::create_directories("results/summary");
    std::filesystem::create_directories("results/chains");

    //Create a file to hold primary and secondary data values
    std::ofstream observables("results/summary/observables_metropolis_chains.dat",std::ios::out);
    observables<<"Temperature MCS AED AMD SH MS"<<std::endl;

    //Initialize lattice and grid points
    model.initialize();

    std::cout<<"Enter length of Monte Carlo Chains: ";
    std::cin>>Nmc;

    if(Nmc<=0){std::cerr<<"Number of Monte Carlo sweeps must be positive."<<std::endl; return;}

    //Define temperature values
    double Temperature[3] = {2.0,2.3,2.6};

    for(double Temp:Temperature){
      //Create lookup table for given temperature
      metropolis.lookup(Temp);

      //Create three unique Monte Carlo chains
      for(int mcgen=1;mcgen<4;mcgen++){
        //Clear everything before use in order to prevent miscalculations
        measurements.energies.clear();
        measurements.magnetizations.clear();

        //Create unique seed
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();;
        generator.seed(seed);

        //Implement spins and calculate total energy and magnetization
        model.discrete_spin_generator(generator,false,"cold");
        model.calculate_total_energy();
        model.calculate_total_magnetization();

        std::cout<<"Temperature = "<<Temp<<"\n";

        //Start measuring time
        auto start = std::chrono::high_resolution_clock::now();

        //Call Metropolis algorithm through Monte Carlo chain
        for(int i=0;i<Nmc;i++){metropolis.sweep(model,generator,measurements);}

        //Calculate average and variance of the equilibrated part
        double en_av = mean_value(measurements.energies,filter);
        double en_var = variance(measurements.energies,en_av,filter);
        double mag_av = mean_value(measurements.magnetizations,filter);
        double mag_var = variance(measurements.magnetizations,mag_av,filter);

        //Calculate average energy and magnetization densities
        summary.energy_density = en_av/model.grid_points;
        summary.magnetization_density = mag_av/model.grid_points;

        //Calculate specific heat and magnetic susceptibility
        calculate_specific_heat(en_var,Temp,model.grid_points,summary);
        calculate_magnetic_susceptibility(mag_var,Temp,model.grid_points,summary);

        std::cout<<"Specific Heat = "<<summary.specific_heat<<"\n"<<"Magnetic Susceptibility = "<<summary.magnetic_susceptibility<<"\n";

        observables<<Temp<<" "<<mcgen<<" "<<summary.energy_density<<" "<<summary.magnetization_density<<" "<<summary.specific_heat<<" "<<summary.magnetic_susceptibility<<std::endl;

        //Export average energy and magnetization densities
        dataexport("results/chains/energies_"+std::to_string(mcgen)+"_"+std::to_string(static_cast<int>(Temp*10)),"dat",measurements.energies,"Times","Energies",filter,100,model.grid_points,true);
        dataexport("results/chains/magnetizations_"+std::to_string(mcgen)+"_"+std::to_string(static_cast<int>(Temp*10)),"dat",measurements.magnetizations,"Times","Magnetizations",filter,100,model.grid_points,true);

        //Stop measuring time and print on the console
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        std::cout<<"Time taken by function: "<<duration.count()<<" microseconds"<<std::endl;
      }
    }
  }
};

int main(){
  MCS obj;
  obj.simulation();
  return 0;
}
