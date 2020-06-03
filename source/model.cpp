#include <fstream>
#include <iostream>
#include <vector>
#include "model.hpp"
#include "grid.hpp"
#include "photons.hpp"
#include "observers.hpp"

Model::Model(Grid *grid, Sources *sources, std::vector<Observer>* observers)
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
    input.ignore(1000000, '=');		input >> R_i;
    input.ignore(1000000, '=');		input >> R_d;
    input.ignore(1000000, '=');		input >> rho_0;
    input.ignore(1000000, '=');		input >> h_0;
    input.ignore(1000000, '=');		input >> R_0;
    input.ignore(1000000, '=');		input >> alpha;
    input.ignore(1000000, '=');		input >> beta;
    // sources parameters
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

    // observers parameters
    double rimage;
    input.ignore(1000000, '=');		input >> rimage;
    std::string typeOfObserversPositions;
    input.ignore(1000000, '=');		input >> typeOfObserversPositions;
    if (typeOfObserversPositions == "MANUAL")
    {
        size_t numberOfObservers;
        input.ignore(1000000, '=');		input >> numberOfObservers;
        observers->reserve(numberOfObservers);
        for (size_t i=0; i!=numberOfObservers; ++i)
        {
            double viewPhi, viewTheta;
            input.ignore(1000000, '=');
            input >> viewPhi >> viewTheta;
            observers->emplace_back(viewPhi*3.1415926/180, viewTheta*3.1415926/180, rimage);
        }
    } else if (typeOfObserversPositions == "PARALLEL") {
        size_t numberOfObservers;
        input.ignore(1000000, '=');		input >> numberOfObservers;
        observers->reserve(numberOfObservers);
        for (size_t i=0; i!=numberOfObservers; ++i)
        {
            double viewTheta;
            input.ignore(1000000, '=');
            input >> viewTheta;
            observers->emplace_back(2*3.141592/numberOfObservers*i, viewTheta*3.1415926/180, rimage);
        }
    } else if (typeOfObserversPositions == "MERIDIAN") {
        size_t numberOfObservers;
        input.ignore(1000000, '=');		input >> numberOfObservers;
        observers->reserve(numberOfObservers);
        for (size_t i=0; i!=numberOfObservers; ++i)
        {
            double viewPhi;
            input.ignore(1000000, '=');
            input >> viewPhi;
            observers->emplace_back(viewPhi*3.1415926/180, 3.141592/(numberOfObservers-1)*i, rimage);
        }
    } else {
        std::cout << "Error in observers positions" << std::endl;
    }

    input.close();
    std::cout << "Parameters\n\nMethod\n fMonteCarlo=" << fMonteCarlo_ << "\n taumin=" << taumin_ << "\n nscat=" << nscat_;
    std::cout << "\n\nMonte Carlo Parameters\n nphotons=" << num_photons_ << "\n iseed=" << iseed_ << "\n\n";
    std::cout << "DGEM Parameters\n PrimaryDirectionsLevel=" << PrimaryDirectionsLevel_ << "\n SecondaryDirectionsLevel=" << SecondaryDirectionsLevel_;
    std::cout << "\n NumOfPrimaryScatterings=" << NumOfPrimaryScatterings_ << "\n NumOfSecondaryScatterings=" << NumOfSecondaryScatterings_;
    std::cout << "\n MonteCarloStart=" << MonteCarloStart_ << "\n\n";
    std::cout << "Physics\n kappa=" << kappa_ << "\n albedo=" << albedo_ << "\n hgg=" << hgg_ << "\n pl=" << pl_ << "\n pc=" << pc_ << "\n sc=" << sc_ << "\n\n";
    std::cout << "Image\n xmax=" << xmax_ << "\n ymax=" << ymax_ << "\n zmax=" << zmax_ << "\n\n";
    std::cout << "Disk\n R_i=" << R_i << "\n R_d=" << R_d << "\n rho_0=" << rho_0 << "\n h_0=" << h_0 << "\n R_0=" << R_0 << "\n alpha=" << alpha;
    std::cout << "\n beta=" << beta;
    // stars
    std::cout << "\nStars\n nstars=" << nstars << "\n";
    for (size_t i=0; i!=nstars; ++i)
        std::cout << " star=" << x[i] << "\t" << y[i] << "\t" << z[i] << "\t" << l[i] << std::endl;

    // observers
    std::cout << "\n\nObservers\n rimage=" << rimage;
    std::cout << "\n nobservers=" << observers->size() << "\n";
    for (size_t i=0; i!=observers->size(); ++i)
        std::cout << " observer=" << (*observers)[i].phi() << "\t" << (*observers)[i].theta() << "\n";

    g2_=hgg_*hgg_;

    grid->Init(*this, R_i, R_d, rho_0, h_0, R_0, alpha, beta, 201, 201, 201);
    sources->Init(nstars, x, y, z, l);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] l;
}
