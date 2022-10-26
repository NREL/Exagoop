#include<iostream>
#include<fstream>
#include "constants.H"
#include "interpolants.H"
#include "hypoplastic.H"
#include "updates.H"

using namespace std;

int main(){
    // input parameters
    GBHypoParameters GBparameters;
    double dt = 1e-6;
    double initial_void_ratio = 0.77;
    double epsilon_dot_11 = -0.1; 
    int output_interval = 1000;
    ofstream datafile;
    datafile.open("data.txt");
    cout<<"Opened datafile for writing, start oedometer test.\n";
    if(datafile.is_open()){
        datafile << "time voidratio sigmaXX sigmaXY sigmaXZ sigmaYY sigmaYZ sigmaZZ epsilonXX \n";
    }
    // initialize oedometer test
    double sigma[NCOMP_TENSOR] = {-10.0, 0.0, 0.0, -0.0, 0.0, -0.0};
    double epsilon[NCOMP_TENSOR] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double omega_dot[NCOMP_TENSOR] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double time = 0;
    double void_ratio = initial_void_ratio;
    int steps = 0;

    // step for loading (compression)
    // stopping criteria is 0.3 MPa for stress XX (T11, sigma_11)
    double epsilon_dot_load[NCOMP_TENSOR] = {epsilon_dot_11, 0.0, 0.0, 0.0, 0.0, 0.0};
    cout<<"start loading\n";
    while(-sigma[XX]<3e5){
        time = time + dt;
        steps ++;
        // update void ratio
        update_void_ratio(void_ratio,epsilon_dot_load,dt);

        // calculate current epsilon: strain
        update_strain(epsilon,epsilon_dot_load,dt);

        // update stress tensor
        GB_hypoplastic(epsilon_dot_load,omega_dot,sigma,dt,void_ratio,GBparameters);

        //output criteria  screen and file
        //cout<< time << " " << void_ratio <<" "<<sigma[XX]<< " "<< sigma[XY]<<" "<<sigma[XZ]<<" "<<sigma[YY]<<" "<<sigma[YZ]<<" "<<sigma[ZZ]<<" "<<epsilon[XX]<<"\n";
        
        
        if(steps % (output_interval/10) == 0){
            cout<< time << " " << void_ratio <<" "<<sigma[XX]<< " "<< sigma[XY]<<" "<<sigma[XZ]<<" "<<sigma[YY]<<" "<<sigma[YZ]<<" "<<sigma[ZZ]<<" "<<epsilon[XX]<<"\n";
        }
        

        if(steps % output_interval == 0){
            datafile<< time << " " << void_ratio <<" "<<sigma[XX]<< " "<< sigma[XY]<<" "<<sigma[XZ]<<" "<<sigma[YY]<<" "<<sigma[YZ]<<" "<<sigma[ZZ]<<" "<<epsilon[XX]<<"\n";
        }

    }

    // step for unloading
    cout<<"start unloading\n";
    double epsilon_dot_unload[NCOMP_TENSOR] = {-epsilon_dot_11, 0.0, 0.0, 0.0, 0.0};
    while(sigma[XX]<-10){
        time = time + dt;
        steps ++;
        // update void ratio
        update_void_ratio(void_ratio,epsilon_dot_unload,dt);

        // calculate current epsilon: strain
        update_strain(epsilon,epsilon_dot_unload,dt);

        // update stress tensor
        GB_hypoplastic(epsilon_dot_unload,omega_dot,sigma,dt,void_ratio,GBparameters);

        //output criteria  screen and file
        if(steps % (output_interval/10) == 0){
            cout<< time << " " << void_ratio <<" "<<sigma[XX]<< " "<< sigma[XY]<<" "<<sigma[XZ]<<" "<<sigma[YY]<<" "<<sigma[YZ]<<" "<<sigma[ZZ]<<" "<<epsilon[XX]<<"\n";
        }

        if(steps % output_interval == 0){
            datafile<< time << " " << void_ratio <<" "<<sigma[XX]<< " "<< sigma[XY]<<" "<<sigma[XZ]<<" "<<sigma[YY]<<" "<<sigma[YZ]<<" "<<sigma[ZZ]<<" "<<epsilon[XX]<<"\n";
        }
    }
    datafile.close();

    cout<<"finished run!\n";
    return 0;


}