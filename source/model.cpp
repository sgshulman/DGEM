#include <iostream>
#include <vector>
#include "model.hpp"
#include "FlaredDisk.hpp"
#include "MathUtils.hpp"
#include "grid.hpp"
#include "SafierWind.hpp"
#include "Sources.hpp"
#include "observers.hpp"
#include "Dust.hpp"
#include "third-party/nlohmann/json.hpp"

namespace
{
    IMatterCPtr parseSafierWind(const nlohmann::json& json)
    {
        return std::make_shared<SafierWind const>(
            json.at("model").get<std::string>().at(0),
            json.at("mOut").get<double>(),
            json.at("mStar").get<double>(),
            json.at("h0").get<double>(),
            json.at("rMin").get<double>(),
            json.at("rMax").get<double>());
    }


    IMatterCPtr parseFlaredDisk(const nlohmann::json& json)
    {
        IMatterCPtr wind;

        if (json.contains("safierWind"))
        {
            wind = parseSafierWind(json.at("safierWind"));
        }

        return std::make_shared<FlaredDisk const>(
            json.at("rinner").get<double>(),
            json.at("router").get<double>(),
            json.at("rho_0").get<double>(),
            json.at("h_0").get<double>(),
            json.at("R_0").get<double>(),
            json.at("alpha").get<double>(),
            json.at("beta").get<double>(),
            wind);
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

    GridCPtr parseGrid(const nlohmann::json& json, double const kappa, IMatterCPtr disk)
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
                    radians(observer.at("phi").get<double>()),
                    radians(observer.at("theta").get<double>()),
                    rimage);
            }
        }

        if (json.contains("parallel"))
        {
            nlohmann::json const& parallellJson = json.at("parallel");
            const auto numberOfObservers = parallellJson.at("numberOfObservers").get<int>();
            const auto viewTheta = parallellJson.at("theta").get<double>();
            observers->reserve(observers->size() + numberOfObservers);

            for (int i=0; i!=numberOfObservers; ++i)
            {
                observers->emplace_back(2*PI/numberOfObservers*i, radians(viewTheta), rimage);
            }
        }

        if (json.contains("meridian"))
        {
            nlohmann::json const& medianJson = json.at("meridian");
            const auto numberOfObservers = medianJson.at("numberOfObservers").get<int>();
            const auto viewPhi = medianJson.at("phi").get<double>();
            observers->reserve(observers->size() + numberOfObservers);

            for (int i=0; i!=numberOfObservers; ++i)
            {
                observers->emplace_back(radians(viewPhi), PI/numberOfObservers*(i + 0.5), rimage);
            }
        }
    }

    SourcesPtr parseSources(nlohmann::json const& json, SourceParameters const& sourceParameters)
    {
        std::vector<PointSource> pointSources;
        pointSources.reserve(json.size());

        for (auto const &star : json)
        {
            pointSources.emplace_back(
                Vector3d{
                    star.at("x").get<double>(),
                    star.at("y").get<double>(),
                    star.at("z").get<double>()},
                star.at("l").get<double>());
        }

        return std::make_shared<Sources>(sourceParameters, std::move(pointSources));
    }
}

#ifdef ENABLE_UNIT_TESTS
#include "tests/ModelTests.inl"
#endif

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
    iseed_ = methodJson.at("iseed").get<int32_t>();

    fMonteCarlo_ = sourceParameters.useMonteCarlo_;
    taumin_ = methodJson.at("taumin").get<double>();
    nscat_ = methodJson.at("nscat").get<uint32_t>();
    SecondaryDirectionsLevel_ = methodJson.at("SecondaryDirectionsLevel").get<uint32_t>();
    NumOfPrimaryScatterings_ = methodJson.at("NumOfPrimaryScatterings").get<uint32_t>();
    NumOfSecondaryScatterings_ = methodJson.at("NumOfSecondaryScatterings").get<uint32_t>();
    MonteCarloStart_ = methodJson.at("MonteCarloStart").get<uint32_t>();

    nlohmann::json const& dustJson = j.at("dust");
    dust_ = parseDust(dustJson);
    IMatterCPtr disk = parseFlaredDisk(j.at("disk"));
    grid_ = parseGrid(j.at("grid"), dustJson.at("kappa").get<double>(), disk);
    sources_ = parseSources(j.at("stars"), sourceParameters);
    parseObservers(observers, j.at("observers"));
}
