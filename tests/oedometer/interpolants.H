#pragma once
#include "constants.H"
#include<iostream>
/**
 * @brief Use for matrix calculation with spin tensor: Get the antisymmetric tensor object based on upper diagonal tensor, 
 * 
 * @param omega_dot : spin tensor
 * @param tens : complete tensor for spin rate
 * @return AMREX_GPU_DEVICE 
 */

void get_antisymmetric_tensor(double omega_dot[NCOMP_TENSOR],
                double tens[AMREX_SPACEDIM*AMREX_SPACEDIM])
{
   int ind=0;
   int start_index=0;
   for(int i=0;i<AMREX_SPACEDIM;i++)
   {
      for(int j=i;j<AMREX_SPACEDIM;j++)
      {
		if(j == i){
			if(omega_dot[start_index+ind]!=0){
				std::cout<<"\nDiagnal element not equal to zero for antisymmetric tensor"<<std::endl;
			}
			else{
				tens[i*AMREX_SPACEDIM+j] = 0;
			}
		}
		else{
			tens[i*AMREX_SPACEDIM+j]=omega_dot[start_index+ind];
        	tens[j*AMREX_SPACEDIM+i]=-omega_dot[start_index+ind];
		}
        ind++;
      }    
   }
}

/**
 * @brief compute inner product of two square matrix with the same order. return (tens1)dot(dens2)
 * 
 * @param tens1 
 * @param tens2 
 * @param product_tensor 
 * @return AMREX_GPU_DEVICE 
 */

void tensor_tensor_dot_product(double tensor1[AMREX_SPACEDIM*AMREX_SPACEDIM],
                double tensor2[AMREX_SPACEDIM*AMREX_SPACEDIM],double product_tensor[AMREX_SPACEDIM*AMREX_SPACEDIM])
{
    for(int i=0; i<AMREX_SPACEDIM; i++){
        for(int j=0; j<AMREX_SPACEDIM; j++){
			product_tensor[i*AMREX_SPACEDIM+j] = 0;
			for(int k=0; k<AMREX_SPACEDIM; k++){
				product_tensor[i*AMREX_SPACEDIM+j] += tensor1[i*AMREX_SPACEDIM+k]*tensor2[j+k*AMREX_SPACEDIM];
			}
        }
    }
}

void get_tensor_from_array(double arr[NCOMP_TENSOR],
                double tens[AMREX_SPACEDIM*AMREX_SPACEDIM])
{
   int ind=0;
   int start_index=0;
   for(int i=0;i<AMREX_SPACEDIM;i++)
   {
      for(int j=i;j<AMREX_SPACEDIM;j++)
      {
        tens[i*AMREX_SPACEDIM+j]=arr[start_index+ind];
        tens[j*AMREX_SPACEDIM+i]=arr[start_index+ind];
        ind++;
      }    
   }
}
