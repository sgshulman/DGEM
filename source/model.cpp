#include <iostream>
#include <vector>
#include "model.hpp"
#include "grid.hpp"
#include "Sources.hpp"
#include "observers.hpp"
#include "Dust.hpp"
#include "third-party/nlohmann/json.hpp"

namespace
{
    FlaredDiskCPtr parseFlaredDisk(const nlohmann::json& json)
    {
        return std::make_shared<FlaredDisk const>(
            json.at("rinner").get<double>(),
            json.at("router").get<double>(),
            json.at("rho_0").get<double>(),
            json.at("h_0").get<double>(),
            json.at("R_0").get<double>(),
            json.at("alpha").get<double>(),
            json.at("beta").get<double>());
    }

    DustCPtr parseDust(const nlohmann::json& json)
    {
        return std::make_shared<Dust>(
            json.at("albedo").get<double>(),
            json.at("hgg").get<double>(),
            json.at("pl").get<double>(),
            json.at("pc").get<double>(),
            json.at("sc").get<double>());
    }

    GridCPtr parseGrid(const nlohmann::json& json, double const kappa, FlaredDiskCPtr disk)
    {
        return std::make_shared<Grid const>(
            json.at("xmax").get<double>(),
            json.at("ymax").get<double>(),
            json.at("zmax").get<double>(),
            kappa,
            201,
            201,
            201,
            std::move(disk));
    }

    void parseObservers(std::vector<Observer>* observers, nlohmann::json const& json)
    {
        auto const rimage = json.at("rimage").get<double>();

        if (json.contains("manual"))
        {
            nlohmann::json const& manualJson = json.at("manual");
            observers->reserve(manualJson.size());
            for (const auto& observer : manualJson)
            {
                observers->emplace_back(
                    observer.at("phi").get<double>()*3.1415926/180,
                    observer.at("theta").get<double>()*3.1415926/180,
                    rimage);
            }
        }

        if (json.contains("parallel"))
        {
            nlohmann::json const& parallellJson = json.at("parallel");
            const auto numberOfObservers = parallellJson.at("numberOfObservers").get<int>();
            const auto viewTheta = parallellJson.at("viewTheta").get<double>();
            observers->reserve(observers->size() + numberOfObservers);

            for (int i=0; i!=numberOfObservers; ++i)
            {
                observers->emplace_back(2*3.141592/numberOfObservers*i, viewTheta*3.1415926/180, rimage);
            }
        }

        if (json.contains("median"))
        {
            nlohmann::json const& medianJson = json.at("median");
            const auto numberOfObservers = medianJson.at("numberOfObservers").get<int>();
            const auto viewPhi = medianJson.at("viewPhi").get<double>();
            observers->reserve(observers->size() + numberOfObservers);

            for (int i=0; i!=numberOfObservers; ++i)
            {
                observers->emplace_back(viewPhi*3.1415926/180, 3.141592/(numberOfObservers-1)*i, rimage);
            }
        }
    }

    SourcesPtr parseSources(nlohmann::json const& json, SourceParameters const& sourceParameters)
    {
        auto *x = new double[json.size()];
        auto *y = new double[json.size()];
        auto *z = new double[json.size()];
        auto *l = new double[json.size()];

            int i = 0;
            for (auto const &star : json)
            {
                x[i] = star.at("x").get<double>();
                y[i] = star.at("y").get<double>();
                z[i] = star.at("z").get<double>();
                l[i] = star.at("l").get<double>();
                ++i;
            }

        auto sources = std::make_shared<Sources>(sourceParameters, i, x, y, z, l);

        delete[] x;
        delete[] y;
        delete[] z;
        delete[] l;

        return sources;
    }
}


Model::Model(std::vector<Observer>* observers)
{
    std::ifstream configurationFile("parameters.json");
    nlohmann::json j;
    configurationFile >> j;

    SourceParameters sourceParameters{};
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
    dust_ = parseDust(dustJson);
    FlaredDiskCPtr disk = parseFlaredDisk(j.at("disk"));
    grid_ = parseGrid(j.at("grid"), dustJson.at("kappa").get<double>(), disk);
    sources_ = parseSources(j.at("stars"), sourceParameters);
    parseObservers(observers, j.at("observers"));
}
