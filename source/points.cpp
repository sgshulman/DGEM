#include <cmath>
#include <iostream>
#include <vector>
#include "model.hpp"
#include "grid.hpp"
#include "observers.hpp"
#include "photons.hpp"
#include "directions.hpp"
#include "Dust.hpp"

Random ran;

void directPhotons(
    SourcesRef sources,
    GridCRef grid,
    std::vector<Observer>* observers,
    unsigned long long const num_photons)
{
    // Direct photon loop.  Loop over sources and weight photons by
    // W=ph*exp(-tau1)/4pi
    for (uint64_t is = 0; is != sources->num(); ++is)
    {
        uint64_t nph = uint64_t(num_photons * (*sources)[is].luminosity() / sources->totlum());

        for (size_t io = 0; io != observers->size(); ++io)
        {
            // Set photon location, grid cell, and direction of observation
            Photon ph((*sources)[is].pos(), Direction3d((*observers)[io].pos()), 1.0, 0);

            // Find optical depth, tau1, to edge of grid along viewing direction
            double tau1 = grid->findOpticalDepth(ph);

            // direct photon weight is exp(-tau1)/4pi
            ph.weight() = nph * exp(-tau1) / 4.0 / PI;

            // bin the photon into the image according to its position and
            // direction of travel.
            (*observers)[io].Bin(ph);
        }
    }
}


int main()
{
    std::vector<Observer> observers;

    // read the model parameters
    Model & model = Model::instance(&observers);
    GridCPtr grid = model.grid();
    SourcesPtr sources = model.sources();

    // Scattered photon loop
    uint64_t totscatt=0;
    uint64_t jcount=0;

    if ( model.fMonteCarlo() )
    {
        // Set up Random
        Random ran( model.iseed() );

        // Loop over sources. nph=number of photons to release from each source
        for (size_t is=0; is!=sources->num(); ++is)
        {
            uint64_t nph=(uint64_t)(model.num_photons() * (*sources)[is].luminosity() / sources->totlum());
            // Loop over nph photons from each source
            for (uint64_t j=0; j<nph; ++j)
            {
                ++jcount;
                if( jcount%10000 == 0)
                {
                    std::cout << jcount << " scattered photons completed" << std::endl;
                }
                // Release photon from point source
                Photon ph((*sources)[is].pos(), 1.0, 1 );

                // Find optical depth, tau1, to edge of grid
                double tau1 = grid->findOpticalDepth( ph );
                if ( tau1 < model.taumin() ) continue;
                double w = 1.0-exp(-tau1);
                ph.weight() = w;

                // Force photon to scatter at optical depth tau before edge of grid
                double tau=-log(1.0-ran.Get()*w );
                // Find scattering location of tau
                grid->movePhotonAtDepth( ph, tau );
                // Photon scatters in grid until it exits (tflag=1) or number
                // of scatterings exceeds a set value (nscatt)
                int tflag = 0;
                while ( !tflag && ( ph.nscat() <= model.nscat() ) )
                {
                    ph.weight() *= model.dust()->albedo();
                    // Do peeling off and project weighted photons into image
                    // учитыается нерассеяшийся свет от каждой точки рассеяния и последующие рассеяния, пока фотон не изыдет
                    for (Observer& observer : observers)
                    {
                        grid->peeloff(ph, observer, model.dust());
                    }

                    // Scatter photon into new direction and update Stokes parameters
                    ph.Stokes( model.dust(), Direction3d(), 0.0, false );
                    ph.nscat()+=1;
                    ++totscatt;
                    if (ph.nscat() > model.nscat()) break;

                    // Find next scattering location
                    tflag = grid->movePhotonAtRandomDepth(ph);
                }
            } // end loop over nph photons
        } // end loop over nsource sources
    } else {
        // Set up directions grid
        Directions pdir( model.PrimaryDirectionsLevel() );
        Directions sdir( model.SecondaryDirectionsLevel() );
        model.num_photons() = pdir.number();
        std::cout << "Directions ready" << std::endl;

        // Loop over sources.
        for (size_t is=0; is!=sources->num(); ++is)
        {
            double w0=(*sources)[is].luminosity()/sources->totlum();
            uint64_t jcount=0;
            // Loop over primary directions
            for (size_t j=0; j!=pdir.number(); ++j)
            {
                ++jcount;
                if( jcount%1000 == 0)
                {
                    std::cout << "Sources: " << is+1 << "/" << sources->num() << ". Directions: " << jcount << "/" << pdir.number() << std::endl;
                }
                // Release photon from point source
                Direction3d direction{ pdir.direction(j) };
                Photon ph0((*sources)[is].pos(), direction, 1.0, 1 );

                // Find optical depth, tau1, to edge of grid
                double tau1 = grid->findOpticalDepth( ph0 );
                if ( tau1 < model.taumin() ) continue;

                double w = (1.0-exp(-tau1)) / model.NumOfPrimaryScatterings() ;
                Vector3d spos = (*sources)[is].pos();
                double tauold = 0.0, tau = 0.0;

                // Loop over scattering dots
                for (size_t s=0; s!=model.NumOfPrimaryScatterings(); ++s)
                {
                    Photon ph( spos, direction, w*w0*pdir.w(j), 1 );
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

                    if (ph.nscat() < model.nscat() ) ph.Scatt( model, sdir, grid, observers );
                }
            }

        }
    }

    directPhotons(sources, grid, &observers, model.num_photons());

    std::cout << "Finishing..." << std::endl;

    // Normalize images by nphotons
    for(size_t cnt=0; cnt!=observers.size(); ++cnt)
    {
        observers[cnt].Normalize( model.num_photons() );
    }

    // put results into output files
    for (size_t cnt = 0; cnt != observers.size(); ++cnt)
    {
        observers[cnt].WriteToMapFiles(true);
    }

    // put general information into file
    std::ofstream observersResultFile("observers.dat");
    for (size_t cnt = 0; cnt != observers.size(); ++cnt)
    {
        observers[cnt].Write(observersResultFile);
    }
    observersResultFile.close();

    return 0;
}
