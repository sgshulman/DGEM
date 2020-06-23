#include <cmath>
#include <iostream>
#include <vector>
#include "model.hpp"
#include "grid.hpp"
#include "observers.hpp"
#include "Photon.hpp"
#include "Directions.hpp"
#include "Dust.hpp"
#include "Sources.hpp"

Random ran;


int main()
{
    std::vector<Observer> observers;

    // read the model parameters
    Model & model = Model::instance(&observers);
    GridCPtr grid = model.grid();
    SourcesPtr sources = model.sources();

    // Scattered photon loop
    if (model.fMonteCarlo())
    {
        for (;;)
        {
            Photon ph{ sources->emitPhoton() };

            if (ph.termination())
            {
                break;
            }

            // Find optical depth, tau1, to edge of grid
            double tau1 = grid->findOpticalDepth(ph);
            if (tau1 < model.taumin())
            {
                continue;
            }
            double w = 1.0 - std::exp(-tau1);
            ph.weight() = w;

            // Force photon to scatter at optical depth tau before edge of grid
            double tau = -std::log(1.0 - ran.Get() * w);
            // Find scattering location of tau
            grid->movePhotonAtDepth(ph, tau);
            // Photon scatters in grid until it exits (tflag=1) or number
            // of scatterings exceeds a set value (nscatt)
            int tflag = 0;
            while (!tflag && (ph.nscat() <= model.nscat()))
            {
                ph.weight() *= model.dust()->albedo();
                // Do peeling off and project weighted photons into image
                // учитыается нерассеяшийся свет от каждой точки рассеяния и последующие рассеяния, пока фотон не изыдет
                for (Observer &observer : observers)
                {
                    grid->peeloff(ph, observer, model.dust());
                }

                // Scatter photon into new direction and update Stokes parameters
                ph.Stokes(model.dust(), Direction3d(), 0.0, false);
                ph.nscat() += 1;
                if (ph.nscat() > model.nscat()) break;

                // Find next scattering location
                tflag = grid->movePhotonAtRandomDepth(ph);
            }
        }
    } else {
        // Set up directions grid
        Directions sdir( model.SecondaryDirectionsLevel() );

        for (;;)
        {
            Photon ph0{ sources->emitPhoton() };

            if (ph0.termination())
            {
                break;
            }

            // Find optical depth, tau1, to edge of grid
            double tau1 = grid->findOpticalDepth( ph0 );
            if ( tau1 < model.taumin() )
            {
                continue;
            }

            double const w = (1.0-exp(-tau1)) / model.NumOfPrimaryScatterings() ;
            double tauold = 0.0, tau = 0.0;
            Vector3d spos = ph0.pos();

            // Loop over scattering dots
            for (size_t s=0; s!=model.NumOfPrimaryScatterings(); ++s)
            {
                Photon ph( spos, ph0.dir(), ph0.weight() * w, 1 );
                // Force photon to scatter at optical depth tau before edge of grid
                tauold = tau;
                tau=-log( 1.0-0.5*w*(2*s+1) );
                // Find scattering location of tau
                grid->movePhotonAtDepth( ph, tau, tauold );
                spos = ph.pos();
                // Photon scatters in grid until it exits (tflag=1) or number
                // of scatterings exceeds a set value (nscatt)

                // Photons scattering
                ph.weight() *= model.dust()->albedo();

                for (Observer& observer : observers)
                {
                    grid->peeloff(ph, observer, model.dust());
                }

                if (ph.nscat() < model.nscat()) ph.Scatt(model, sdir, grid, observers);
            }
        }
    }

    sources->directPhotons(grid, &observers);

    std::cout << "Finishing..." << std::endl;

    // Normalize images
    for(size_t cnt=0; cnt!=observers.size(); ++cnt)
    {
        observers[cnt].normalize(sources->num_photons());
    }

    // put results into output files
    for (size_t cnt = 0; cnt != observers.size(); ++cnt)
    {
        observers[cnt].writeToMapFiles(true);
    }

    // put general information into file
    std::ofstream observersResultFile("observers.dat");
    for (size_t cnt = 0; cnt != observers.size(); ++cnt)
    {
        observers[cnt].write(observersResultFile);
    }
    observersResultFile.close();

    return 0;
}
