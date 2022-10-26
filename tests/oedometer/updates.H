#pragma once
#include "constants.H"
void update_strain(double strain[NCOMP_TENSOR], double strain_rate[NCOMP_TENSOR], double dt){
    for(int i = 0; i < NCOMP_TENSOR; i++){
        strain[i] += strain_rate[i]*dt;
    }
}

void update_void_ratio(double& void_ratio, double strain_rate[NCOMP_TENSOR], double dt){
    void_ratio += (strain_rate[XX]+strain_rate[YY]+strain_rate[ZZ])* (1 +void_ratio)* dt;
}