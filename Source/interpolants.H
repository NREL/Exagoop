#ifndef PARTICLEBASIS_H_
#define PARTICLEBASIS_H_

#include<constants.H>
#include <AMReX.H>
#include<mpm_particle_container.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real spherical_gaussian(amrex::Real xi[3],amrex::Real xp[3],amrex::Real r0)
{
  amrex::Real f;
  amrex::Real r2;

  r2=(xi[0]-xp[0])*(xi[0]-xp[0])+(xi[1]-xp[1])*(xi[1]-xp[1])+(xi[2]-xp[2])*(xi[2]-xp[2]);
  f=std::exp(-r2/(r0*r0))*std::pow(PI,-1.5)*std::pow(r0,-3.0);

  return(f);
}

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
                amrex::Array4<amrex::Real> nodaldata,
				int comp)
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
        amrex::Real cubicspline_1d_der(int shapefunctiontype,amrex:: Real zi)
{
	amrex::Real value=0.0;
	if(shapefunctiontype==1)	//Boundary node
	{
		if(zi>=-2 && zi<=-1)
		{
			value=zi*zi/2.0+2*zi+2.0;
		}
		else if(zi>=-1 && zi<=0)
		{
			value=-zi*zi/2.0+1.0;
		}
		else if(zi>=0 && zi<=1)
		{
			value=zi*zi/2.0-1.0;
		}
		else if(zi>=1 && zi<=2.0)
		{
			value=-zi*zi/2.0+2*zi-2.0;
		}
	}
	else if(shapefunctiontype==2)	//Near Boundary node
	{
		if(zi>=-1 && zi<=0)
		{
			value=-zi*zi-2.0*zi;
		}
		else if(zi>=0 && zi<=1)
		{
			value=3.0*zi*zi/2.0-2.0*zi;
		}
		else if(zi>=1 && zi<=2.0)
		{
			value=-zi*zi/2.0+2*zi-2.0;
		}
	}
	else if(shapefunctiontype==3)	//Interior node
	{
		if(zi>=-2 && zi<=-1)
		{
			value=zi*zi/2.0+2*zi+2.0;
		}
		else if(zi>=-1 && zi<=0)
		{
			value=-3.0*zi*zi/2.0-2.0*zi;
		}
		else if(zi>=0 && zi<=1)
		{
			value=3.0*zi*zi/2.0-2.0*zi;
		}
		else if(zi>=1 && zi<=2.0)
		{
			value=-zi*zi/2.0+2*zi-2.0;
		}
	}
	else if(shapefunctiontype==4)	//Interior node
	{
		if(zi>=-2 && zi<=-1)
		{
			value=zi*zi/2.0+2*zi+2.0;
		}
		else if(zi>=-1 && zi<=0)
		{
			value=-3.0*zi*zi/2.0-2.0*zi;
		}
		else if(zi>=0 && zi<=1)
		{
			value=zi*zi-2.0*zi;
		}
	}
	else
	{
#ifndef AMREX_USE_GPU
		amrex::Print()<<"\n Shapefunction = "<<shapefunctiontype;
#endif
		amrex::Abort("\n Incorrect shape function type in cubicspline_1d_der");
	}
	return value;

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real cubicspline_1d(int shapefunctiontype,amrex:: Real zi)
{
	amrex::Real value=0.0;
	if(shapefunctiontype==1)	//Boundary node
	{
		if(zi>=-2 && zi<=-1)
		{
			value=zi*zi*zi/6.0+zi*zi+2.0*zi+4.0/3.0;
		}
		else if(zi>=-1 && zi<=0)
		{
			value=-zi*zi*zi/6.0+zi+1.0;
		}
		else if(zi>=0 && zi<=1)
		{
			value=zi*zi*zi/6.0-zi+1.0;
		}
		else if(zi>=1 && zi<=2.0)
		{
			value=-zi*zi*zi/6.0+zi*zi-2.0*zi+4.0/3.0;
		}
	}
	else if(shapefunctiontype==2)	//Near Boundary node
	{
		if(zi>=-1 && zi<=0)
		{
			value=-zi*zi*zi/3.0-zi*zi+2.0/3.0;
		}
		else if(zi>=0 && zi<=1)
		{
			value=zi*zi*zi/2.0-zi*zi+2.0/3.0;
		}
		else if(zi>=1 && zi<=2.0)
		{
			value=-zi*zi*zi/6.0+zi*zi-2.0*zi+4.0/3.0;
		}
	}
	else if(shapefunctiontype==3)	//Interior node
	{
		if(zi>=-2 && zi<=-1)
		{
			value=zi*zi*zi/6.0+zi*zi+2.0*zi+4.0/3.0;
		}
		else if(zi>=-1 && zi<=0)
		{
			value=-zi*zi*zi/2.0-zi*zi+2.0/3.0;
		}
		else if(zi>=0 && zi<=1)
		{
			value=zi*zi*zi/2.0-zi*zi+2.0/3.0;
		}
		else if(zi>=1 && zi<=2.0)
		{
			value=-zi*zi*zi/6.0+zi*zi-2.0*zi+4.0/3.0;
		}
	}
	else if(shapefunctiontype==4)	//Interior node
	{
		if(zi>=-2 && zi<=-1)
		{
			value=zi*zi*zi/6.0+zi*zi+2.0*zi+4.0/3.0;
		}
		else if(zi>=-1 && zi<=0)
		{
			value=-zi*zi*zi/2.0-zi*zi+2.0/3.0;
		}
		else if(zi>=0 && zi<=1)
		{
			value=zi*zi*zi/3.0-zi*zi+2.0/3.0;
		}
		else
                {
#ifndef AMREX_USE_GPU
			amrex::Print()<<"\n Error in shape func 4: "<<zi;
#endif
		}
	}
	else
	{
#ifndef AMREX_USE_GPU
		amrex::Print()<<"\n Shapefunction = "<<shapefunctiontype;
#endif
		amrex::Abort("\n Incorrect shape function type in cubicspline_1d");
	}
	return value;

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real cubic_interp(amrex::Real xp[AMREX_SPACEDIM],
        		int i, int j, int k,
        		int lmin,int mmin,int nmin,int lmax,int mmax, int nmax,
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo,
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
                amrex::Array4<amrex::Real> nodaldata,
				int comp,const int* lo,const int* hi)
{

    amrex::Real rx,ry,rz;
    amrex::Real lval,mval,nval;
    amrex::Real value=zero;
    int shapetypex,shapetypey,shapetypez;

    for(int n=nmin;n<nmax;n++)
    {
    	for(int m=mmin;m<mmax;m++)
    	{
    		for(int l=lmin;l<lmax;l++)
    		{
    			if(i+l==lo[0]) shapetypex=1;
    			else if(i+l==lo[0]+1) shapetypex=2;
    			else if(i+l==hi[0]+1) shapetypex=1;
    			else if(i+l==hi[0]) shapetypex=4;
    			else shapetypex=3;

    			if(j+m==lo[1]) shapetypey=1;
    			else if(j+m==lo[1]+1) shapetypey=2;
    			else if(j+m==hi[1]+1) shapetypey=1;
    			else if(j+m==hi[1]) shapetypey=4;
    			else shapetypey=3;

    			if(k+n==lo[2]) shapetypez=1;
    			else if(k+n==lo[2]+1) shapetypez=2;
    			else if(k+n==hi[2]+1) shapetypez=1;
    			else if(k+n==hi[2]) shapetypez=4;
    			else shapetypez=3;

    			rx = (xp[XDIR]-(plo[XDIR]+(i+l)*dx[XDIR]))/dx[XDIR];
    			ry  = (xp[YDIR]-(plo[YDIR]+(j+m)*dx[YDIR]))/dx[YDIR];
    			rz  = (xp[ZDIR]-(plo[ZDIR]+(k+n)*dx[ZDIR]))/dx[ZDIR];

    			lval=cubicspline_1d(shapetypex,rx);
    			mval=cubicspline_1d(shapetypey,ry);
    			nval=cubicspline_1d(shapetypez,rz);
    			value += lval*mval*nval* nodaldata(i+l,j+m,k+n,comp);

    		}
    	}
    }

    return(value);
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real basisval(int l,int m,int n, int i, int j, int k,
               amrex::Real xp[AMREX_SPACEDIM], 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo, 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
				int order_scheme,
				const int *lo,
				const int *hi)
{
    amrex::Real lval,mval,nval;
    amrex::Real rx,ry,rz;

    lval=0.0;
    mval=0.0;
    nval=0.0;
    
    int shapefunctype;

    rx 	= (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    ry  = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    rz  = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];

    if(order_scheme==1)
    {
    	lval=(l==0)?(one-rx):rx;
    	mval=(m==0)?(one-ry):ry;
    	nval=(n==0)?(one-rz):rz;
    }
    else if(order_scheme==3)
    {
    	//Checking for i node
    	if(i==lo[0])
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==0)?1:((l==1)?2:((l==2)?3:-1));
    		rx=(l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==lo[0]+1)
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?1:((l==0)?2:((l==1)?3:(l==2)?3:-1));
    		if(shapefunctype==-1)
    		{
    			amrex::Abort("\n Incorrect shape function at i ==lo[0]+1");
    		}
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==lo[0]+2)
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?2:((l==0)?3:((l==1)?3:(l==2)?3:-1));
    		if(shapefunctype==-1)
    		{
    			amrex::Abort("\n Incorrect shape function at i ==lo[0]+1");
    		}
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==hi[0])
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?3:((l==0)?4:((l==1)?1:-1));
    		if(shapefunctype==-1)
    		{
    			amrex::Abort("\n Incorrect shape function at i ==hi[0]-1");
    		}
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):-1000));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==hi[0]-1)
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?3:(((l==0)?3:((l==1)?4:(l==2)?1:-1)));
    		if(shapefunctype==-1)
    		    		{
    		    			amrex::Abort("\n Incorrect shape function at i ==hi[0]-1");
    		    		}
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-1-(1-rx):-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==hi[0]-2)
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?3:((l==0)?3:((l==1)?3:((l==2)?4:-1)));
    		if(shapefunctype==-1)
    		    		{
    		    			amrex::Abort("\n Incorrect shape function at i ==hi[0]-2");
    		    		}
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=3;
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}

    	//Checking for j node

    	if(j==lo[1])
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==0)?1:((m==1)?2:((m==2)?3:-1));
    		ry=(m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==lo[1]+1)
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?1:((m==0)?2:((m==1)?3:(m==2)?3:-1));
    		if(shapefunctype==-1)
    		{
    			amrex::Abort("\n Incorrect shape function at i ==lo[1]+1");
    		}
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==lo[1]+2)
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?2:((m==0)?3:((m==1)?3:(m==2)?3:-1));
    		if(shapefunctype==-1)
    		{
    			amrex::Abort("\n Incorrect shape function at i ==lo[1]+2");
    		}
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==hi[1])
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?3:((m==0)?4:((m==1)?1:-1));
    		if(shapefunctype==-1)
    		{
    			amrex::Abort("\n Incorrect shape function at i ==hi[0]");
    		}
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):-1000));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==hi[1]-1)
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?3:(((m==0)?3:((m==1)?4:(m==2)?1:-1)));
    		if(shapefunctype==-1)
    		{
    			amrex::Abort("\n Incorrect shape function at i ==hi[1]-1");
    		}
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-1-(1-ry):-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==hi[1]-2)
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?3:((m==0)?3:((m==1)?3:((m==2)?4:-1)));
    		if(shapefunctype==-1)
    		{
    			amrex::Abort("\n Incorrect shape function at i ==hi[0]-2");
    		}
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=3;
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}

    	//Checking for k node

		if(k==lo[2])
		{
			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
			shapefunctype=(n==0)?1:((n==1)?2:((n==2)?3:-1));
			rz=(n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000));
			nval=cubicspline_1d(shapefunctype,rz);
		}
		else if(k==lo[2]+1)
		{
			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
			shapefunctype=(n==-1)?1:((n==0)?2:((n==1)?3:(n==2)?3:-1));
			if(shapefunctype==-1)
			{
				amrex::Abort("\n Incorrect shape function at i ==lo[1]+1");
			}
			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
			nval=cubicspline_1d(shapefunctype,rz);
		}
		else if(k==lo[2]+2)
		{
			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
			shapefunctype=(n==-1)?2:((n==0)?3:((n==1)?3:(n==2)?3:-1));
			if(shapefunctype==-1)
			{
				amrex::Abort("\n Incorrect shape function at i ==lo[1]+2");
			}
			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
			nval=cubicspline_1d(shapefunctype,rz);
		}
		else if(k==hi[2])
		{
			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
			shapefunctype=(n==-1)?3:((n==0)?4:((n==1)?1:-1));
			if(shapefunctype==-1)
			{
				amrex::Abort("\n Incorrect shape function at i ==hi[0]");
			}
			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):-1000));
			nval=cubicspline_1d(shapefunctype,rz);
		}
		else if(k==hi[2]-1)
		{
			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
			shapefunctype=(n==-1)?3:(((n==0)?3:((n==1)?4:(n==2)?1:-1)));
			if(shapefunctype==-1)
			{
				amrex::Abort("\n Incorrect shape function at i ==hi[1]-1");
			}
			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-1-(1-rz):-1000)));
			nval=cubicspline_1d(shapefunctype,rz);
		}
		else if(k==hi[2]-2)
		{
			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
			shapefunctype=(n==-1)?3:((n==0)?3:((n==1)?3:((n==2)?4:-1)));
			if(shapefunctype==-1)
			{
				amrex::Abort("\n Incorrect shape function at i ==hi[0]-2");
			}
			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
			nval=cubicspline_1d(shapefunctype,rz);
		}
		else
		{
			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
			shapefunctype=3;
			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
			nval=cubicspline_1d(shapefunctype,rz);
		}

    }
    else
    {
        amrex::Abort("order scheme not implemented yet\n");
    }
    return(lval*mval*nval);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real basisvalder(int dir, int l,int m,int n, int i, int j, int k, 
                amrex::Real xp[AMREX_SPACEDIM], 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo, 
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
				int order_scheme,const int *lo,
				const int *hi)
{
    amrex::Real lval,mval,nval;
    amrex::Real rx,ry,rz;
    amrex::Real dxinv[AMREX_SPACEDIM];

    int shapefunctype;

    dxinv[XDIR]=one/dx[XDIR];
    dxinv[YDIR]=one/dx[YDIR];
    dxinv[ZDIR]=one/dx[ZDIR];

    rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    ry  = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    rz  = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];

    if(order_scheme==1)
    {
    	lval=(l==0)?(one-rx):rx;
    	mval=(m==0)?(one-ry):ry;
    	nval=(n==0)?(one-rz):rz;

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
    	//amrex::Print()<<"\nEntered basisvalder"<<lval<<" "<<mval<<" "<<nval;
    }
    else if(order_scheme==3)
    {
    	//Checking for i node
    	if(i==lo[0])
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==0)?1:((l==1)?2:((l==2)?3:-1));
    		rx=(l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==lo[0]+1)
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?1:((l==0)?2:((l==1)?3:(l==2)?3:-1));
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==lo[0]+2)
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?2:((l==0)?3:((l==1)?3:(l==2)?3:-1));
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==hi[0])
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?3:((l==0)?4:((l==1)?1:-1));
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):-1000));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==hi[0]-1)
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?3:((l==0)?3:((l==1)?4:((l==2)?1:-1)));
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else if(i==hi[0]-2)
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=(l==-1)?3:((l==0)?3:((l==1)?3:((l==2)?4:-1)));
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}
    	else
    	{
    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		shapefunctype=3;
    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		lval=cubicspline_1d(shapefunctype,rx);
    	}



    	//Checking for j node

    	if(j==lo[1])
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==0)?1:((m==1)?2:((m==2)?3:-1));
    		ry=(m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==lo[1]+1)
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?1:((m==0)?2:((m==1)?3:(m==2)?3:-1));
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==lo[1]+2)
		{
			ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
			shapefunctype=(m==-1)?2:((m==0)?3:((m==1)?3:(m==2)?3:-1));
			ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
			mval=cubicspline_1d(shapefunctype,ry);
		}
    	else if(j==hi[1])
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?3:((m==0)?4:((m==1)?1:-1));
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):-1000));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==hi[1]-1)
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?3:((m==0)?3:((m==1)?4:((m==2)?1:-1)));
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else if(j==hi[1]-2)
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=(m==-1)?3:((m==0)?3:((m==1)?3:((m==2)?4:-1)));
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}
    	else
    	{
    		ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    		shapefunctype=3;
    		ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    		mval=cubicspline_1d(shapefunctype,ry);
    	}

    	//Checking for k node

    	if(k==lo[2])
    	{
    		rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    		shapefunctype=(n==0)?1:((n==1)?2:((n==2)?3:-1));
    		rz=(n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000));
    		nval=cubicspline_1d(shapefunctype,rz);
    	}
    	else if(k==lo[2]+1)
    	{
    		rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    		shapefunctype=(n==-1)?1:((n==0)?2:((n==1)?3:(n==2)?3:-1));
    		rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    		nval=cubicspline_1d(shapefunctype,rz);
    	}
    	else if(k==lo[2]+2)
    	{
    		rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    		shapefunctype=(n==-1)?2:((n==0)?3:((n==1)?3:(n==2)?3:-1));
    		rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    		nval=cubicspline_1d(shapefunctype,rz);
    	}
    	else if(k==hi[2])
    	{
    		rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    		shapefunctype=(n==-1)?3:((n==0)?4:((n==1)?1:-1));
    		rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):-1000));
    		nval=cubicspline_1d(shapefunctype,rz);
    	}
    	else if(k==hi[2]-1)
    	{
    		rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    		shapefunctype=(n==-1)?3:((n==0)?3:((n==1)?4:((n==2)?1:-1)));
    		rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    		nval=cubicspline_1d(shapefunctype,rz);
    	}
    	else if(k==hi[2]-2)
    	{
    		rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    		shapefunctype=(n==-1)?3:((n==0)?3:((n==1)?3:((n==2)?4:-1)));
    		rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    		nval=cubicspline_1d(shapefunctype,rz);
    	}
    	else
    	{
    		rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    		shapefunctype=3;
    		rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    		nval=cubicspline_1d(shapefunctype,rz);
    	}


    	if(dir==0)
    	{
    		if(i==lo[0])
    		{
    			rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    			shapefunctype=(l==0)?1:((l==1)?2:((l==2)?3:-1));
    			rx=(l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000));
    			lval=cubicspline_1d_der(shapefunctype,rx)*dxinv[XDIR];
    		}
    		else if(i==lo[0]+1)
    		{
    			rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    			shapefunctype=(l==-1)?1:((l==0)?2:((l==1)?3:(l==2)?3:-1));
    			rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    			lval=cubicspline_1d_der(shapefunctype,rx)*dxinv[XDIR];
    		}
    		else if(i==lo[0]+2)
			{
				rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
				shapefunctype=(l==-1)?2:((l==0)?3:((l==1)?3:(l==2)?3:-1));
				rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
				lval=cubicspline_1d_der(shapefunctype,rx)*dxinv[XDIR];
			}
    		else if(i==hi[0])
    		{
    			rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    			shapefunctype=(l==-1)?3:((l==0)?4:((l==1)?1:-1));
    			rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):-1000));
    			lval=cubicspline_1d_der(shapefunctype,rx)*dxinv[XDIR];
    		}
    		else if(i==hi[0]-1)
    		{
    			rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    			shapefunctype=(l==-1)?3:((l==0)?3:((l==1)?4:((l==2)?1:-1)));
    			rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    			lval=cubicspline_1d_der(shapefunctype,rx)*dxinv[XDIR];
    		}
    		else if(i==hi[0]-2)
    		    	{
    		    		rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    		    		shapefunctype=(l==-1)?3:((l==0)?3:((l==1)?3:((l==2)?4:-1)));
    		    		rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    		    		lval=cubicspline_1d_der(shapefunctype,rx)*dxinv[XDIR];
    		    	}
    		else
    		{
    			rx = (xp[XDIR]-(plo[XDIR]+i*dx[XDIR]))/dx[XDIR];
    			shapefunctype=3;
    			rx=(l==-1)?1+rx:((l==0)?rx:((l==1)?-(1-rx):((l==2)?-(1-rx)-1:-1000)));
    			lval=cubicspline_1d_der(shapefunctype,rx)*dxinv[XDIR];
    		}
    	}
    	else if(dir==1)
    	{
    		if(j==lo[1])
    		{
    			ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    			shapefunctype=(m==0)?1:((m==1)?2:((m==2)?3:-1));
    			ry=(m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000));
    			mval=cubicspline_1d_der(shapefunctype,ry)*dxinv[YDIR];
    		}
    		else if(j==lo[1]+1)
    		{
    			ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    			shapefunctype=(m==-1)?1:((m==0)?2:((m==1)?3:(m==2)?3:-1));
    			ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    			mval=cubicspline_1d_der(shapefunctype,ry)*dxinv[YDIR];
    		}
    		else if(j==lo[1]+2)
    		{
    			ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    			shapefunctype=(m==-1)?2:((m==0)?3:((m==1)?3:(m==2)?3:-1));
    			ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    			mval=cubicspline_1d_der(shapefunctype,ry)*dxinv[YDIR];
    		}
    		else if(j==hi[1])
    		{
    			ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    			shapefunctype=(m==-1)?3:((m==0)?4:((m==1)?1:-1));
    			ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):-1000));
    			mval=cubicspline_1d_der(shapefunctype,ry)*dxinv[YDIR];

    		}
    		else if(j==hi[1]-1)
    		{
    			ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    			shapefunctype=(m==-1)?3:((m==0)?3:((m==1)?4:((m==2)?1:-1)));
    			ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    			mval=cubicspline_1d_der(shapefunctype,ry)*dxinv[YDIR];
    		}
    		else if(j==hi[1]-2)
    		{
    			ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    			shapefunctype=(m==-1)?3:((m==0)?3:((m==1)?3:((m==2)?4:-1)));
    			ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    			mval=cubicspline_1d_der(shapefunctype,ry)*dxinv[YDIR];
    		}
    		else
    		{
    			ry = (xp[YDIR]-(plo[YDIR]+j*dx[YDIR]))/dx[YDIR];
    			shapefunctype=3;
    			ry=(m==-1)?1+ry:((m==0)?ry:((m==1)?-(1-ry):((m==2)?-(1-ry)-1:-1000)));
    			mval=cubicspline_1d_der(shapefunctype,ry)*dxinv[YDIR];
    		}
    	}
    	else if(dir==2)
    	{
    		if(k==lo[2])
    		{
    			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    			shapefunctype=(n==0)?1:((n==1)?2:((n==2)?3:-1));
    			rz=(n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000));
    			nval=cubicspline_1d_der(shapefunctype,rz)*dxinv[ZDIR];
    		}
    		else if(k==lo[2]+1)
    		{
    			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    			shapefunctype=(n==-1)?1:((n==0)?2:((n==1)?3:(n==2)?3:-1));
    			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    			nval=cubicspline_1d_der(shapefunctype,rz)*dxinv[ZDIR];
    		}
    		else if(k==lo[2]+2)
    		{
    			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    			shapefunctype=(n==-1)?2:((n==0)?3:((n==1)?3:(n==2)?3:-1));
    			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    			nval=cubicspline_1d_der(shapefunctype,rz)*dxinv[ZDIR];
    		}
    		else if(k==hi[2])
    		{
    			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    			shapefunctype=(n==-1)?3:((n==0)?4:((n==1)?1:-1));
    			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):-1000));
    			nval=cubicspline_1d_der(shapefunctype,rz)*dxinv[ZDIR];
    		}
    		else if(k==hi[2]-1)
    		{
    			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    			shapefunctype=(n==-1)?3:((n==0)?3:((n==1)?4:((n==2)?1:-1)));
    			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    			nval=cubicspline_1d_der(shapefunctype,rz)*dxinv[ZDIR];
    		}
    		else if(k==hi[2]-2)
    		{
    			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    			shapefunctype=(n==-1)?3:((n==0)?3:((n==1)?3:((n==2)?4:-1)));
    			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    			nval=cubicspline_1d_der(shapefunctype,rz)*dxinv[ZDIR];
    		}
    		else
    		{
    			rz = (xp[ZDIR]-(plo[ZDIR]+k*dx[ZDIR]))/dx[ZDIR];
    			shapefunctype=3;
    			rz=(n==-1)?1+rz:((n==0)?rz:((n==1)?-(1-rz):((n==2)?-(1-rz)-1:-1000)));
    			nval=cubicspline_1d_der(shapefunctype,rz)*dxinv[ZDIR];
    		}
    	}
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void get_tensor_from_array(amrex::Real arr[NCOMP_TENSOR],
                amrex::Real tens[AMREX_SPACEDIM*AMREX_SPACEDIM])
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

/**
 * @brief Use for matrix calculation with spin tensor: Get the antisymmetric tensor object based on upper diagonal tensor, 
 * 
 * @param omega_dot : spin tensor
 * @param tens : complete tensor for spin rate
 * @return AMREX_GPU_DEVICE 
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void get_antisymmetric_tensor(amrex::Real omega_dot[NCOMP_TENSOR],
                amrex::Real tens[AMREX_SPACEDIM*AMREX_SPACEDIM])
{
   int ind=0;
   int start_index=0;
   for(int i=0;i<AMREX_SPACEDIM;i++)
   {
      for(int j=i;j<AMREX_SPACEDIM;j++)
      {
		if(j == i){
			if(omega_dot[start_index+ind]!=0){
				amrex::Abort("\nDiagnal element not equal to zero for antisymmetric tensor");
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void tensor_tensor_dot_product(amrex::Real tensor1[AMREX_SPACEDIM*AMREX_SPACEDIM],
                amrex::Real tensor2[AMREX_SPACEDIM*AMREX_SPACEDIM],amrex::Real product_tensor[AMREX_SPACEDIM*AMREX_SPACEDIM])
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



#endif
