#ifndef _INOUTPUT_HPP_
#define _INOUTPUT_HPP_

#include <fstream>
#include <iostream>
#include "model.hpp"
#include "photons.hpp"

// pictures 
class PICTURES
{
	public:
		PICTURES(void) : rimage_(0), Nx_(0), Ny_(0)
		{
			f_ = nullptr;
			q_ = nullptr;
			u_ = nullptr;
		}
		void Init(double rimage,uint32_t Nx=200, uint32_t Ny=200) 
		{
			rimage_=rimage;
			Nx_=Nx; 
			Ny_=Ny;
			f_ = new double[ Nx*Ny ];
			q_ = new double[ Nx*Ny ];
			u_ = new double[ Nx*Ny ];

			for (size_t cnt=0; cnt!=Nx*Ny; ++cnt)
			{
				f_[cnt]=0.0;
				q_[cnt]=0.0;
				u_[cnt]=0.0;
			}
		}
		~PICTURES()
		{
			delete[] f_;
			delete[] q_;
			delete[] u_;
		}
		// place phton on the images
		void Bin( PHOTON ph )
		{
			DIRECTION o=ph.dir();
        	double yimage=rimage_+ph.z()*o.sint()-ph.y()*o.cost()*o.sinp()-ph.x()*o.cost()*o.cosp();
		    double ximage=rimage_+ph.y()*o.cosp()-ph.x()*o.sinp();
		    int64_t xl=int(Nx_*ximage/(2.0*rimage_));
		    int64_t yl=int(Ny_*yimage/(2.0*rimage_));
			// place weighted photon into image location
        	if((xl>=0) && (yl >= 0) && ((uint32_t)xl < Nx_) && ((uint32_t)yl < Ny_) ) 
        		f_[xl+yl*Nx_] += ph.weight();     
        }
        void Bin( PHOTON ph, double f, double q, double u )
		{
			DIRECTION o=ph.dir();
        	double yimage=rimage_+ph.z()*o.sint()-ph.y()*o.cost()*o.sinp()-ph.x()*o.cost()*o.cosp();
		    double ximage=rimage_+ph.y()*o.cosp()-ph.x()*o.sinp();
		    int64_t xl=int(Nx_*ximage/(2.0*rimage_));
		    int64_t yl=int(Ny_*yimage/(2.0*rimage_));
			// place weighted photon into image location
        	if((xl>=0) && (yl >= 0) && ((uint32_t)xl < Nx_) && ((uint32_t)yl < Ny_) ) 
        	{ 
        		f_[xl+yl*Nx_] += f;
        		q_[xl+yl*Nx_] += q;
        		u_[xl+yl*Nx_] += u;
        	}     
        }
            
	void Norm( int num )
	{
		double imtot = 0.0;
		for (size_t i=0; i!=Nx_*Ny_; ++i)
		{
	       f_[i]=f_[i]/num;
	       q_[i]=q_[i]/num;
	       u_[i]=u_[i]/num;
	       imtot += f_[i];
		}
    	std::cout << "4*pi*imtot = " << imtot*4*PI << std::endl;
    }	
    void Write( int key ) const
    {
		char fname[20];
		sprintf(fname, "fimage%2.2i.dat", key);
		std::ofstream f(fname);
        for (size_t x=0; x!=Nx_; ++x)
		{
			for (size_t y=0; y!=Ny_; ++y)
				f << f_[y+x*Ny_] << "\t";
			f << "\n";
		}
        f.close();
		sprintf(fname, "qimage%2.2i.dat", key);
		std::ofstream q(fname);
        for (size_t x=0; x!=Nx_; ++x)
		{
			for (size_t y=0; y!=Ny_; ++y)
				q << q_[y+x*Ny_] << "\t";
			q << "\n";
		}
        q.close();
          	
        sprintf(fname, "uimage%2.2i.dat", key);
		std::ofstream u(fname);
        for (size_t x=0; x!=Nx_; ++x)
		{
			for (size_t y=0; y!=Ny_; ++y)
				u << u_[y+x*Ny_] << "\t";
			u << "\n";
		}
        u.close();  
    }
    void Sum()
    {
		double fsum=0, qsum=0, usum=0;
		for (size_t x=0; x!=Nx_; ++x)
		{
			for (size_t y=0; y!=Ny_; ++y)
			{
				fsum += f_[y+x*Ny_];
				qsum += q_[y+x*Ny_];
				usum += u_[y+x*Ny_];
			}
		}		
		std::cout << "F = " << fsum << "\nQ = " << qsum << "\nU = " << usum << "\n\n"; 
	}
	
	private:
		double rimage_;
		uint32_t Nx_, Ny_;
		double *f_;
		double *q_;
		double *u_;

		PICTURES ( PICTURES const & pics );
		PICTURES & operator =( PICTURES const & pics );
};

#endif
