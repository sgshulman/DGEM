#include <iostream>
#include <vector>
#include "model.hpp"
#include "grid.hpp"
#include "Sources.hpp"
#include "observers.hpp"
#include "Dust.hpp"
#include "third-party/nlohmann/json.hpp"

Model::Model(std::vector<Observer>* observers)
{
    std::ifstream configurationFile("parameters.json");
    nlohmann::json j;
    configurationFile >> j;

    SourceParameters sourceParameters;
    nlohmann::json const& methodJson = j.at("method parameters");
    sourceParameters.useMonteCarlo_ = methodJson.at("fMonteCarlo").get<bool>();
    sourceParameters.num_photons_ = methodJson.at("nphotons").get<uint64_t>();
    sourceParameters.PrimaryDirectionsLevel_ = methodJson.at("PrimaryDirectionsLevel").get<uint32_t>();
    sourceParameters.seed_ = methodJson.at("iseed").get<int32_t>();

    fMonteCarlo_ = sourceParameters.useMonteCarlo_;
    taumin_ = methodJson.at("taumin").get<double>();
    nscat_ = methodJson.at("nscat").get<uint32_t>();
    SecondaryDirectionsLevel_ = methodJson.at("SecondaryDirectionsLevel").get<uint32_t>();
    NumOfPrimaryScatterings_ = methodJson.at("NumOfPrimaryScatterings").get<uint32_t>();
    NumOfSecondaryScatterings_ = methodJson.at("NumOfSecondaryScatterings").get<uint32_t>();
    MonteCarloStart_ = methodJson.at("MonteCarloStart").get<uint32_t>();

    nlohmann::json const& dustJson = j.at("dust");

    dust_ = std::make_shared<Dust>(
        dustJson.at("albedo").get<double>(),
        dustJson.at("hgg").get<double>(),
        dustJson.at("pl").get<double>(),
        dustJson.at("pc").get<double>(),
        dustJson.at("sc").get<double>());

    nlohmann::json const& diskJson = j.at("disk");

    FlaredDiskCPtr disk = std::make_shared<FlaredDisk const>(
        diskJson.at("rinner").get<double>(),
        diskJson.at("router").get<double>(),
        diskJson.at("rho_0").get<double>(),
        diskJson.at("h_0").get<double>(),
        diskJson.at("R_0").get<double>(),
        diskJson.at("alpha").get<double>(),
        diskJson.at("beta").get<double>());

    nlohmann::json const& gridJson = j.at("grid");

    grid_ = std::make_shared<Grid const>(
        gridJson.at("xmax").get<double>(),
        gridJson.at("ymax").get<double>(),
        gridJson.at("zmax").get<double>(),
        dustJson.at("kappa").get<double>(),
        201,
        201,
        201,
        disk);

    nlohmann::json const& starsJson = j.at("stars");

    auto *x = new double[starsJson.size()];
    auto *y = new double[starsJson.size()];
    auto *z = new double[starsJson.size()];
    auto *l = new double[starsJson.size()];

    {
        int i = 0;
        for (auto const &star : starsJson)
        {
            x[i] = star.at("x").get<double>();
            y[i] = star.at("y").get<double>();
            z[i] = star.at("z").get<double>();
            l[i] = star.at("l").get<double>();
            ++i;
        }

        sources_ = std::make_shared<Sources>(sourceParameters, i, x, y, z, l);
    }

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] l;

    nlohmann::json const& observersJson = j.at("observers");
    auto const rimage = observersJson.at("rimage").get<double>();

    if (observersJson.contains("manual"))
    {
        nlohmann::json const& manualJson = observersJson.at("manual");
        observers->reserve(manualJson.size());
        for (const auto& observer : manualJson)
        {
            observers->emplace_back(
                observer.at("phi").get<double>()*3.1415926/180,
                observer.at("theta").get<double>()*3.1415926/180,
                rimage);
        }
    }

    if (observersJson.contains("parallel"))
    {
        nlohmann::json const& parallellJson = observersJson.at("parallel");
        const auto numberOfObservers = parallellJson.at("numberOfObservers").get<int>();
        const auto viewTheta = parallellJson.at("viewTheta").get<double>();
        observers->reserve(observers->size() + numberOfObservers);

        for (int i=0; i!=numberOfObservers; ++i)
        {
            observers->emplace_back(2*3.141592/numberOfObservers*i, viewTheta*3.1415926/180, rimage);
        }
    }

    if (observersJson.contains("median"))
    {
        nlohmann::json const& medianJson = observersJson.at("median");
        const auto numberOfObservers = medianJson.at("numberOfObservers").get<int>();
        const auto viewPhi = medianJson.at("viewPhi").get<double>();
        observers->reserve(observers->size() + numberOfObservers);

        for (int i=0; i!=numberOfObservers; ++i)
        {
            observers->emplace_back(viewPhi*3.1415926/180, 3.141592/(numberOfObservers-1)*i, rimage);
        }
    }
}
