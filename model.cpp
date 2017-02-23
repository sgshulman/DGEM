#include <fstream>
#include <iostream>
#include "model.hpp"
#include "grid.hpp"
#include "photons.hpp"

MODEL::MODEL (GRID *grid, SOURCES *sources) 
{
	// disk parameters
	double R_i, R_d;
	double rho_0;
	double h_0, R_0;
	double alpha, beta;
        
    // sources parameters
    uint32_t nstars;
    double *x, *y, *z, *l;
    std::ifstream input("params.par");
	input.ignore(1000000, '=');		input >> fMonteCarlo_;
	input.ignore(1000000, '=');		input >> fTauHash_;
	input.ignore(1000000, '=');		input >> taumin_;
	input.ignore(1000000, '=');		input >> nscat_;
	input.ignore(1000000, '=');		input >> num_photons_;
	input.ignore(1000000, '=');		input >> iseed_;
	input.ignore(1000000, '=');		input >> PrimaryDirectionsLevel_;
	input.ignore(1000000, '=');		input >> SecondaryDirectionsLevel_;
	input.ignore(1000000, '=');		input >> NumOfPrimaryScatterings_;
	input.ignore(1000000, '=');		input >> NumOfSecondaryScatterings_;
	input.ignore(1000000, '=');		input >> MonteCarloStart_;
	input.ignore(1000000, '=');		input >> kappa_;
	input.ignore(1000000, '=');		input >> albedo_;
	input.ignore(1000000, '=');		input >> hgg_;
	input.ignore(1000000, '=');		input >> pl_;
	input.ignore(1000000, '=');		input >> pc_;
	input.ignore(1000000, '=');		input >> sc_;
	input.ignore(1000000, '=');		input >> xmax_;
	input.ignore(1000000, '=');		input >> ymax_;
	input.ignore(1000000, '=');		input >> zmax_;
	input.ignore(1000000, '=');		input >> rimage_;
	input.ignore(1000000, '=');		input >> viewthet_;
	input.ignore(1000000, '=');		input >> viewphi_;
	input.ignore(1000000, '=');		input >> R_i;
	input.ignore(1000000, '=');		input >> R_d;
	input.ignore(1000000, '=');		input >> rho_0;
	input.ignore(1000000, '=');		input >> h_0;
	input.ignore(1000000, '=');		input >> R_0;
	input.ignore(1000000, '=');		input >> alpha;
	input.ignore(1000000, '=');		input >> beta;
	input.ignore(1000000, '=');		input >> nstars;
	x = new double[nstars];
	y = new double[nstars];
	z = new double[nstars];
	l = new double[nstars];
	for (size_t i=0; i!=nstars; ++i)
	{
		input.ignore(1000000, '=');
		input >> x[i] >> y[i] >> z[i] >> l[i];
	}
	input.close();
	std::cout << "Parameters\n\nMethod\n fMonteCarlo=" << fMonteCarlo_ << "\n fTauHash=" << fTauHash_ << "\n taumin=" << taumin_ << "\n nscat=" << nscat_;
	std::cout << "\n\nMonte Carlo Parameters\n nphotons=" << num_photons_ << "\n iseed=" << iseed_ << "\n\n";
	std::cout << "DGEM Parameters\n PrimaryDirectionsLevel=" << PrimaryDirectionsLevel_ << "\n SecondaryDirectionsLevel=" << SecondaryDirectionsLevel_;
	std::cout << "\n NumOfPrimaryScatterings=" << NumOfPrimaryScatterings_ << "\n NumOfSecondaryScatterings=" << NumOfSecondaryScatterings_;
	std::cout << "\n MonteCarloStart=" << MonteCarloStart_ << "\n\n";
	std::cout << "Physics\n kappa=" << kappa_ << "\n albedo=" << albedo_ << "\n hgg=" << hgg_ << "\n pl=" << pl_ << "\n pc=" << pc_ << "\n sc=" << sc_ << "\n\n";
	std::cout << "Image\n xmax=" << xmax_ << "\n ymax=" << ymax_ << "\n zmax=" << zmax_ << "\n rimage=" << rimage_ << "\n\n";
	std::cout << "Observer position\n viewthet=" << viewthet_ << "\n viewphi=" << viewphi_ << "\n\n";
    std::cout << "Disk\n R_i=" << R_i << "\n R_d=" << R_d << "\n rho_0=" << rho_0 << "\n h_0=" << h_0 << "\n R_0=" << R_0 << "\n alpha=" << alpha;
    std::cout << "\n beta=" << beta;
    std::cout << "\nStars\n nstars=" << nstars << "\n";
    for (size_t i=0; i!=nstars; ++i)
	    std::cout << " star=" << x[i] << "\t" << y[i] << "\t" << z[i] << "\t" << l[i] << std::endl;
			    
	g2_=hgg_*hgg_;
			
	grid->Init(*this, R_i, R_d, rho_0, h_0, R_0, alpha, beta, 201, 201, 201);
	sources->Init(nstars, x, y, z, l);
			
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] l;
}
