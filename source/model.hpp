#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <cstdint>
#include <vector>

#include "Predefines.hpp"

double const PI	= 3.141592653589793238;

// Random number generator 
// L’Ecuyer, P. 1988, Communications of the ACM, vol. 31, pp. 742,774.
class Random
{
    public:
		explicit Random (int32_t iseed=-1556) : iseed_(iseed) { } ;
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


// model parameters
class Model
{
	public:
		static Model& instance (std::vector<Observer> *observers)
		{
			static Model mod(observers) ;
			return mod ;
		}
		bool fMonteCarlo() const
		{	return fMonteCarlo_;	}

		double taumin() const
		{	return taumin_;	}
		uint32_t nscat() const
		{	return nscat_;	}
		int32_t iseed() const
		{	return iseed_;	}
		uint32_t SecondaryDirectionsLevel() const
		{	return SecondaryDirectionsLevel_;	}
		uint32_t MonteCarloStart() const
		{	return MonteCarloStart_;	}
		uint32_t NumOfPrimaryScatterings() const
		{	return NumOfPrimaryScatterings_;	}
		uint32_t NumOfSecondaryScatterings() const
		{	return NumOfSecondaryScatterings_;	}
		GridCPtr grid() const
        {   return grid_; }
		DustCPtr dust() const
        {   return dust_; }
        SourcesPtr sources()
        {   return sources_; }
	private:
		bool fMonteCarlo_;
		double taumin_;
		int32_t iseed_;

        uint32_t nscat_;
		uint32_t MonteCarloStart_;
		uint32_t SecondaryDirectionsLevel_;
        uint32_t NumOfPrimaryScatterings_;
		uint32_t NumOfSecondaryScatterings_;

		DustCPtr dust_;
        GridCPtr grid_;
        SourcesPtr sources_;

		Model(std::vector<Observer> *observers);
		Model(Model const &);
		Model& operator=( Model const &);
};

#endif
