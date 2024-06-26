#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include "model.hpp"
#include "CartesianGrid.hpp"
#include "IDust.hpp"
#include "IRandomGenerator.hpp"
#include "Observer.hpp"
#include "Photon.hpp"
#include "Directions.hpp"
#include "Sources.hpp"

namespace
{

void processMonteCarloPhotons(
    Model const& model,
    std::vector<Observer>& observers,
    IGridCRef grid,
    SourcesRef sources,
    IRandomGenerator *ran)
{
    for (;;)
    {
        double randomValue{ 0.0 };
        Photon ph{ sources->emitRandomPhoton(grid, ran, &randomValue) };

        if (ph.termination())
        {
            return;
        }

        // Find optical depth, tau1, to edge of grid
        double tau1 = grid->findOpticalDepth(ph);
        if (tau1 < model.taumin())
        {
            ran->Skip();
            continue;
        }
        double w = 1.0 - std::exp(-tau1);
        ph.weight() = w;

        // Force photon to scatter at optical depth tau before edge of grid
        double tau = -std::log(1.0 - randomValue * w);
        // Find scattering location of tau
        grid->movePhotonAtDepth(ph, tau, 0.0);

        // Photon scatters inside grid until leaves it or the number
        // of scatterings exceeds a set value (nscatt)
        bool insideGrid = true;
        while (insideGrid)
        {
            ph.weight() *= model.dust()->albedo();
            // Do peeling off and project weighted photons into image
            for (Observer &observer : observers)
            {
                grid->peeloff(ph, observer, model.dust());
            }

            if (ph.nscat() >= model.nscat())
            {
                break;
            }

            // Scatter photon into new direction and update Stokes parameters
            ph.Stokes(model.dust(), Direction3d(), 0.0, false, ran);
            ph.nscat() += 1;

            // Find next scattering location
            insideGrid = grid->movePhotonAtRandomDepth(ph, ran);
        }
        ran->Skip();
    }
}


int run(const std::string& parametersFileName)
{
    auto const initTime = std::chrono::high_resolution_clock::now();

    std::vector<Observer> observers;

    // read the model parameters
    Model& model = Model::instance(&observers, parametersFileName);
    IGridCPtr grid = model.grid();
    SourcesPtr sources = model.sources();

    IRandomGenerator* ran{ model.createRandomGenerator() };
    if (ran != nullptr)
    {
        std::cout << "Random generator type:\t" << ran->GetConfiguration() << std::endl;
    }

    std::cout << "Matter mass:\t" << grid->computeMatterMass() << "\t Solar Masses" << std::endl;
    sources->writeObserversOpticalDepths(grid, &observers);

    auto startTime = std::chrono::high_resolution_clock::now();

    // Scattered photon loop
    if (model.fMonteCarlo())
    {
        processMonteCarloPhotons(model, observers, grid, sources, ran);
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

                for (std::uint32_t s = 0; s != model.NumOfPrimaryScatterings(); ++s)
                {
                    Photon ph(spos, sCellId, ph0.dir(), ph0.weight() * w, 1);
                    // Force photon to scatter at optical depth tau before edge of grid
                    tauold = tau;
                    tau = -log(1.0 - 0.5 * w * (2 * s + 1));

                    // Find scattering location of tau
                    grid->movePhotonAtDepth(ph, tau, tauold);
                    spos = ph.pos();
                    sCellId = ph.cellId();

                    // Photons scattering
                    ph.weight() *= model.dust()->albedo();

                    for (Observer &observer : observers)
                    {
                        grid->peeloff(ph, observer, model.dust());
                    }

                    if (ph.nscat() < model.nscat())
                    {
                        ph.Scatt(model, sdir, grid, observers, ran);
                    }

                    if (ran != nullptr)
                    {
                        ran->Skip();
                    }
                }
            }
        }
        else
        {
            double const sqrtPiN = std::sqrt(PI / static_cast<double>(sources->num_photons()));
            double const base = 2. * sqrtPiN / (1 - sqrtPiN);
            // pessimistic estimation of scattering number, as we use it only for minimal optical depth estimation
            double const nScatteringsRev = std::log(1 + base) / std::log(std::sqrt(3.) * grid->max());
            double const minWeight = 1. - std::exp(-model.taumin() * nScatteringsRev);
            std::unique_ptr<IRandomGenerator> dgemRandom{ model.createDgemRandomGenerator() };

            for (;;)
            {
                Photon ph0{sources->emitDgemPhoton(grid)};

                if (ph0.termination())
                {
                    break;
                }

                double tau1 = grid->findOpticalDepth(ph0);
                if (tau1 < model.taumin())
                {
                    continue;
                }

                Vector3d spos = ph0.pos();

                // skip empty inner regions
                grid->movePhotonAtDepth(ph0, std::numeric_limits<double>::epsilon(), 0.0);
                double oldR = std::max(model.defaultStarRadius(), (spos - ph0.pos()).norm());
                double baseMultiplier = dgemRandom ? dgemRandom->Get() : 1.;

                while (grid->inside(ph0) && ph0.weight() > minWeight)
                {
                    Photon ph{ ph0 };

                    // estimate tau
                    double const r = oldR * (1. + base * baseMultiplier);
                    baseMultiplier = 1.;
                    Vector3d oldpos = ph0.pos();
                    double const tau = grid->movePhotonAtDistance(ph0, r - oldR);
                    ph0.weight() *= std::exp(-tau);

                    // Photons scattering
                    if (tau > model.taumin() * nScatteringsRev * nScatteringsRev)
                    {
                        // Find scattering location of tau
                        double const tauScattering = -log(0.5 + 0.5 * exp(-tau));
                        grid->movePhotonAtDepth(ph, tauScattering, 0.0);

                        ph.weight() *= model.dust()->albedo() * (1.0 - std::exp(-tau));

                        switch (model.dgemBinType())
                        {
                            case (DgemBinType::POINT):
                                for (Observer &observer : observers)
                                {
                                     grid->peeloff(ph, observer, model.dust());
                                }
                                break;
                            case (DgemBinType::LINE):
                                for (Observer &observer : observers)
                                {
                                     grid->peeloff(ph, observer, model.dust(), oldpos, ph0.pos());
                                }
                                break;
                            case (DgemBinType::HEX_LINES):
                                for (Observer &observer : observers)
                                {
                                     grid->peeloffHex(ph, observer, model.dust(), oldpos, ph0.pos());
                                }
                                break;
                            default:
                                ;
                        }

                        if (ph.nscat() < model.nscat())
                        {
                            ph.Scatt(model, sdir, grid, observers, ran);
                        }
                    }

                    if (ran != nullptr)
                    {
                        ran->Skip();
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
        observers[cnt].writeToMapFiles(model.writeScatterings());
    }

    // put general information into file
    std::ofstream observersResultFile("observers.dat");
    observersResultFile.precision(14);
    for (std::uint64_t cnt = 0; cnt != observers.size(); ++cnt)
    {
        observers[cnt].write(observersResultFile);
    }
    observersResultFile.close();

    if (ran != nullptr)
    {
        ran->save();
        delete ran;
    }

    auto const endTime = std::chrono::high_resolution_clock::now();
    std::cout.precision(4);
    std::cout << "Initialization time:\t" << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(startTime - initTime).count()) * 0.001 << "s\n"
              << "Computation time:\t" << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()) * 0.001 << "s" << std::endl;

    {
        std::ofstream timefile("time.txt");
        timefile.precision(4);

        timefile << "Initialization time:    "
                 << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(startTime - initTime).count()) * 0.001 << "s\n"
                 << "Computation time:       "
                 << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()) * 0.001 << "s"
                 << std::endl;

        timefile.close();
    }

    return 0;
}
} // namespace

int main(int argc, char *argv[])
{
    try
    {
        if (argc > 2)
        {
            std::cerr << argv[0] << " may use only 1 argument, but " << argc-1 << " were provided.\n"
                << "Usage:\n" << argv[0] << " [-v | --version | PARAMETERS_FILE]" << std::endl;
        }
        else if (argc == 2)
        {
            std::string parameter = argv[1];
            if (parameter == "-v" || parameter == "--version")
            {
                std::cout << "Version 0.5\n"
                    << "    Implement Quasi Monte-Carlo method" << std::endl;
            }
            else
            {
                return run(parameter);
            }
        }
        else
        {
            return run("parameters.json");
        }
    }
    catch (std::exception const& ex) {
        std::cerr << "Terminated after throwing an exception:\n" << ex.what() << std::endl;
    }
    return 1;
}
