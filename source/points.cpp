#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include "model.hpp"
#include "CartesianGrid.hpp"
#include "Observer.hpp"
#include "Photon.hpp"
#include "Random.hpp"
#include "Directions.hpp"
#include "IDust.hpp"
#include "Sources.hpp"


int run(int argc, char *argv[])
{
    auto const initTime = std::chrono::high_resolution_clock::now();

    std::vector<Observer> observers;

    // read the model parameters
    Model& model = Model::instance(&observers, argc == 2 ? argv[1] : "parameters.json");
    IGridCPtr grid = model.grid();
    SourcesPtr sources = model.sources();
    Random ran{ model.createRandomGenerator() };

    std::cout << "Matter mass:\t" << grid->computeMatterMass() << "\t Solar Masses" << std::endl;
    sources->writeObserversOpticalDepths(grid, &observers);

    auto startTime = std::chrono::high_resolution_clock::now();

    // Scattered photon loop
    if (model.fMonteCarlo())
    {
        for (;;)
        {
            Photon ph{ sources->emitRandomPhoton(grid, &ran) };

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
            grid->movePhotonAtDepth(ph, tau, 0.0);
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
                ph.Stokes(model.dust(), Direction3d(), 0.0, false, &ran);
                ph.nscat() += 1;
                if (ph.nscat() > model.nscat()) break;

                // Find next scattering location
                tflag = grid->movePhotonAtRandomDepth(ph, &ran);
            }
        }
    } else {
        // Set up directions grid
        Directions sdir( model.SecondaryDirectionsLevel(), model.useHEALPixGrid() );
        startTime = std::chrono::high_resolution_clock::now();

        // Loop over scattering dots
        if (model.NumOfPrimaryScatterings() > 0)
        {
            for (;;)
            {
                Photon ph0{sources->emitDgemPhoton(grid)};

                if (ph0.termination())
                {
                    break;
                }

                // Find optical depth, tau1, to edge of grid
                double tau1 = grid->findOpticalDepth(ph0);
                if (tau1 < model.taumin())
                {
                    continue;
                }

                double const w = (1.0 - exp(-tau1)) / model.NumOfPrimaryScatterings();

                double tauold = 0.0, tau = 0.0;
                Vector3d spos = ph0.pos();
                std::uint64_t sCellId = ph0.cellId();

                for (std::uint64_t s = 0; s != model.NumOfPrimaryScatterings(); ++s)
                {
                    Photon ph(spos, sCellId, ph0.dir(), ph0.weight() * w, 1);
                    // Force photon to scatter at optical depth tau before edge of grid
                    tauold = tau;
                    tau = -log(1.0 - 0.5 * w * (2 * s + 1));

                    // Find scattering location of tau
                    grid->movePhotonAtDepth(ph, tau, tauold);
                    spos = ph.pos();
                    sCellId = ph.cellId();

                    // Photon scatters in grid until it exits (tflag=1) or number
                    // of scatterings exceeds a set value (nscatt)

                    // Photons scattering
                    ph.weight() *= model.dust()->albedo();

                    for (Observer &observer : observers)
                    {
                        grid->peeloff(ph, observer, model.dust());
                    }

                    if (ph.nscat() < model.nscat()) ph.Scatt(model, sdir, grid, observers, &ran);
                }
            }
        }
        else
        {
            double const sqrtPiN = std::sqrt(PI / sources->num_photons());
            double const base = 1. + 2. * sqrtPiN / (1 - sqrtPiN);
            // pessimistic estimation of scattering number, as we use it only for minimal optical depth estimation
            double const nScatteringsRev = std::log(base) / std::log(std::sqrt(3.) * grid->max());
            double const minWeight = 1. - std::exp(-model.taumin() * nScatteringsRev);
            std::cout << minWeight << std::endl;

            for (;;)
            {
                Photon ph0{sources->emitDgemPhoton(grid)};

                if (ph0.termination())
                {
                    break;
                }

                Vector3d spos = ph0.pos();

                // skip empty inner regions
                grid->movePhotonAtDepth(ph0, std::numeric_limits<double>::epsilon(), 0.0);
                double oldR = std::max(model.defaultStarRadius(), (spos - ph0.pos()).norm());

                while (grid->inside(ph0) && ph0.weight() > minWeight)
                {
                    Photon ph{ ph0 };

                    // estimate tau
                    double const r = oldR * base;
                    Vector3d oldpos = ph0.pos();
                    double const tau = grid->movePhotonAtDistance(ph0, r - oldR);
                    ph0.weight() *= std::exp(-tau);

                    // Photons scattering
                    if (tau > model.taumin() * nScatteringsRev)
                    {
                        // Find scattering location of tau
                        double const tauScattering = -log(0.5 + 0.5 * exp(-tau));
                        grid->movePhotonAtDepth(ph, tauScattering, 0.0);

                        ph.weight() *= model.dust()->albedo() * (1.0 - std::exp(-tau));

                        for (Observer &observer : observers)
                        {
                            grid->peeloff(ph, observer, model.dust(), oldpos, ph0.pos());
                        }

                        if (ph.nscat() < model.nscat()) ph.Scatt(model, sdir, grid, observers, &ran);
                    }
                    oldR = r;
                }
            }
        }
    }

    sources->directPhotons(grid, &observers);

    std::cout << "Finishing..." << std::endl;

    // Normalize images
    for(std::uint64_t cnt=0; cnt!=observers.size(); ++cnt)
    {
        observers[cnt].normalize(sources->num_photons());
    }

    // put results into output files
    for (std::uint64_t cnt = 0; cnt != observers.size(); ++cnt)
    {
        observers[cnt].writeToMapFiles(model.writeScatterings(), model.nscat());
    }

    // put general information into file
    std::ofstream observersResultFile("observers.dat");
    for (std::uint64_t cnt = 0; cnt != observers.size(); ++cnt)
    {
        observers[cnt].write(observersResultFile);
    }
    observersResultFile.close();

    ran.save();

    auto const endTime = std::chrono::high_resolution_clock::now();
    std::cout.precision(4);
    std::cout << "Initialization time:\t" << std::chrono::duration_cast<std::chrono::milliseconds>(startTime - initTime).count() * 0.001 << "s\n"
              << "Computation time:\t" << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() * 0.001 << "s" << std::endl;

    {
        std::ofstream timefile("time.txt");
        timefile.precision(4);

        timefile << "Initialization time:    "
                 << std::chrono::duration_cast<std::chrono::milliseconds>(startTime - initTime).count() * 0.001 << "s\n"
                 << "Computation time:       "
                 << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() * 0.001 << "s"
                 << std::endl;

        timefile.close();
    }

    return 0;
}


int main(int argc, char *argv[])
{
    try
    {
        return run(argc, argv);
    }
    catch (std::exception const& ex) {
        std::cerr << "Terminated after throwing an exception:\n" << ex.what() << std::endl;
    }
    return 0;
}
