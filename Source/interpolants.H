#ifndef PARTICLEBASIS_H_
#define PARTICLEBASIS_H_

#include<constants.H>
#include <AMReX.H>
#include<mpm_particle_container.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real hat1d(amrex::Real xi,amrex::Real xp,amrex::Real hatsize)
{
    amrex::Real zeta=xi-xp;
    amrex::Real halfhatsize=hatsize*half;
    amrex::Real funcval=zero;
    amrex::Real maxval=one/halfhatsize;
    amrex::Real slope=four/(hatsize*hatsize);

    if(amrex::Math::abs(zeta) < halfhatsize)
    {
        if(zeta>0.0)
        {
           funcval = -slope*zeta + maxval;       
        }
        else
        {
           funcval =  slope*zeta + maxval;
        }
    }
    return(funcval);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real hat3d(amrex::Real xi[AMREX_SPACEDIM], amrex::Real xp[AMREX_SPACEDIM],
                amrex::Real hatsize[AMREX_SPACEDIM])
{
    amrex::Real funcval=one;
    funcval *= hat1d(xi[XDIR],xp[XDIR],hatsize[XDIR]);
    funcval *= hat1d(xi[YDIR],xp[YDIR],hatsize[YDIR]);
    funcval *= hat1d(xi[ZDIR],xp[ZDIR],hatsize[ZDIR]);

    return(funcval);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real bilin_interp(amrex::Real xp[AMREX_SPACEDIM], int i, int j, int k, 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo, 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
                amrex::Array4<amrex::Real> nodaldata, int comp) 
{

    amrex::Real alpha,beta,gama;
    amrex::Real value=zero;

    alpha = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    beta  = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    gama  = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];

    value += (1.0-alpha) * (1-beta) * (1-gama) * nodaldata(i,j,k,comp);
    value += alpha       * (1-beta) * (1-gama) * nodaldata(i+1,j,k,comp);
    value += (1.0-alpha) * beta     * (1-gama) * nodaldata(i,j+1,k,comp);
    value += alpha       * beta     * (1-gama) * nodaldata(i+1,j+1,k,comp);
    
    value += (1.0-alpha) * (1-beta)  * gama  * nodaldata(i,j,k+1,comp);
    value += alpha       * (1-beta)  * gama  * nodaldata(i+1,j,k+1,comp);
    value += (1.0-alpha) * beta      * gama  * nodaldata(i,j+1,k+1,comp);
    value += alpha       * beta      * gama  * nodaldata(i+1,j+1,k+1,comp);
    
    return(value);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real basisval(int l,int m,int n, int i, int j, int k,
               amrex::Real xp[AMREX_SPACEDIM], 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo, 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    amrex::Real lval,mval,nval;
    amrex::Real alpha,beta,gama;
    
    alpha = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    beta  = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    gama  = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];

    lval=(l==0)?(one-alpha):alpha;
    mval=(m==0)?(one-beta):beta;
    nval=(n==0)?(one-gama):gama;

    return(lval*mval*nval);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real basisvalder(int dir, int l,int m,int n, int i, int j, int k, 
                amrex::Real xp[AMREX_SPACEDIM], 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo, 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    amrex::Real lval,mval,nval;
    amrex::Real alpha,beta,gama;
    amrex::Real dxinv[AMREX_SPACEDIM];

    dxinv[XDIR]=one/dx[XDIR];
    dxinv[YDIR]=one/dx[YDIR];
    dxinv[ZDIR]=one/dx[ZDIR];
    
    alpha = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    beta  = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    gama  = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];

    lval=(l==0)?(one-alpha):alpha;
    mval=(m==0)?(one-beta):beta;
    nval=(n==0)?(one-gama):gama;

    if(dir==0)
    {
        lval=(l==0)?-dxinv[XDIR]:dxinv[XDIR];
    }
    else if(dir==1)
    {
        mval=(m==0)?-dxinv[YDIR]:dxinv[YDIR];
    }
    else
    {
        nval=(n==0)?-dxinv[ZDIR]:dxinv[ZDIR];
    }

    return(lval*mval*nval);
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void get_tensor(MPMParticleContainer::ParticleType &p,int start_index,
                amrex::Real tens[AMREX_SPACEDIM*AMREX_SPACEDIM])
{
   int ind=0;
   for(int i=0;i<AMREX_SPACEDIM;i++)
   {
      for(int j=i;j<AMREX_SPACEDIM;j++)
      {
        tens[i*AMREX_SPACEDIM+j]=p.rdata(start_index+ind);
        tens[j*AMREX_SPACEDIM+i]=p.rdata(start_index+ind);
        ind++;
      }    
   }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void tensor_vector_pdt(amrex::Real tens[AMREX_SPACEDIM*AMREX_SPACEDIM],
                amrex::Real vect[AMREX_SPACEDIM],amrex::Real tensvect[AMREX_SPACEDIM])
{
    for(int i=0;i<AMREX_SPACEDIM;i++)
    {
        tensvect[i]=zero;
        for(int j=0;j<AMREX_SPACEDIM;j++)
        {
            tensvect[i] += tens[i*AMREX_SPACEDIM+j]*vect[j];
        }
    }
}



#endif