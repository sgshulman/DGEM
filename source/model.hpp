#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <cstdint>
#include <vector>

double const PI	= 3.1415926;

class Grid ;
class Pictures ;
class Observer;
class Photon ;
class Source ;
class Sources ;
class Directions ;  // directions grid

// Random number generator 
// L’Ecuyer, P. 1988, Communications of the ACM, vol. 31, pp. 742,774.
class Random
{
    public:
		Random (int32_t iseed=-1556) : iseed_(iseed) { } ;
		double Get()
		{
			  int32_t IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV,j,k;
			  double ran2,AM,EPS,RNMX;
			  IM1=2147483563;IM2=2147483399;AM=1./IM1;IMM1=IM1-1;
			  IA1=40014;IA2=40692;IQ1=53668;IQ2=52774;IR1=12211;IR2=3791;
			  NTAB=32;NDIV=1+IMM1/NTAB;EPS=1.2e-7;RNMX=1.-EPS;
			  static int32_t idum2=123456789,iy=0;
			  static int32_t iv[32]={0};
			  
			  if (iseed_ <= 0) 
			  {
				iseed_=(-iseed_ > 1) ? -iseed_ : 1;
				idum2=iseed_;
				for (j=NTAB+7; j>=0;--j)
				{
				  k=iseed_/IQ1;
				  iseed_=IA1*(iseed_-k*IQ1)-k*IR1;
				  if (iseed_ < 0) iseed_=iseed_+IM1;
				  if (j < NTAB) iv[j]=iseed_;
				}
				iy=iv[0];
			  }
			  k=iseed_/IQ1;
			  iseed_=IA1*(iseed_-k*IQ1)-k*IR1;
			  if (iseed_ < 0) iseed_=iseed_+IM1;
			  k=idum2/IQ2;
			  idum2=IA2*(idum2-k*IQ2)-k*IR2;
			  if (idum2 < 0) idum2=idum2+IM2;
			  j=iy/NDIV;
			  iy=iv[j]-idum2;
			  iv[j]=iseed_;
			  if(iy < 1)iy=iy+IMM1;
			  ran2=(AM*iy < RNMX) ? AM*iy : RNMX;
			  return (double)ran2;
		}
	private:
		int iseed_;
};
extern Random ran;

// model parameters
class Model
{
	public:
		static Model & instance (Grid *grid, Sources *sources, std::vector<Observer> *observers)
		{
			static Model mod(grid, sources, observers) ;
			return mod ;
		}
		bool fMonteCarlo() const
		{	return fMonteCarlo_;	}
		double taumin() const
		{	return taumin_;	}
		uint32_t nscat() const
		{	return nscat_;	}
        uint64_t& num_photons()
		{	return num_photons_;	}
        uint64_t num_photons() const
		{	return num_photons_;	}
		int32_t iseed() const
		{	return iseed_;	}
		uint32_t PrimaryDirectionsLevel() const
		{	return PrimaryDirectionsLevel_;	}
		uint32_t SecondaryDirectionsLevel() const
		{	return SecondaryDirectionsLevel_;	}
		uint32_t MonteCarloStart() const
		{	return MonteCarloStart_;	}
		uint32_t NumOfPrimaryScatterings() const
		{	return NumOfPrimaryScatterings_;	}
		uint32_t NumOfSecondaryScatterings() const
		{	return NumOfSecondaryScatterings_;	}
		double kappa() const
		{	return kappa_;	}
		double albedo() const
		{	return albedo_;	}
		double hgg() const
		{	return hgg_;	}
		double g2() const
		{	return g2_;		}
		double pl() const
		{	return pl_;		}
		double pc() const
		{	return pc_;		}
		double sc() const
		{	return sc_;		}
		double xmax() const
		{	return xmax_;	}
		double ymax() const
		{	return ymax_;	}
		double zmax() const
		{	return zmax_;	}
		double rimage() const
		{	return rimage_;	}

	private:
		bool fMonteCarlo_;
		double taumin_;
		uint32_t nscat_;
		uint64_t num_photons_;
		int32_t iseed_;
		uint32_t PrimaryDirectionsLevel_;
		uint32_t SecondaryDirectionsLevel_;
		uint32_t MonteCarloStart_;
		uint32_t NumOfPrimaryScatterings_, NumOfSecondaryScatterings_;	
		double kappa_;
		double albedo_;
		double hgg_, g2_;
		double pl_;
		double pc_;
		double sc_;
		double xmax_, ymax_, zmax_;
		double rimage_;

		Model (Grid *grid, Sources *sources, std::vector<Observer> *observers);
		Model (Model const &);
		Model & operator =( Model const &);
};

#endif
