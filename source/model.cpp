#include <iostream>
#include <vector>
#include "model.hpp"
#include "DebugUtils.hpp"
#include "Dust.hpp"
#include "FlaredDisk.hpp"
#include "MathUtils.hpp"
#include "MatterArray.hpp"
#include "MatterTranslation.hpp"
#include "grid.hpp"
#include "SafierWind.hpp"
#include "Sources.hpp"
#include "SphereEnvelope.hpp"
#include "observers.hpp"
#include "third-party/nlohmann/json.hpp"

namespace
{
    MatterTranslationCPtr parseTranslation(const nlohmann::json& json)
    {
        return std::make_shared<MatterTranslation const>(
            json.contains("precession") ? radians(json.at("precession").get<double>()) : 0.0,
            json.contains("nutation") ? radians(json.at("nutation").get<double>()) : 0.0,
            json.contains("intrinsicRotation") ? radians(json.at("intrinsicRotation").get<double>()) : 0.0,
            Vector3d{
                json.contains("x") ? json.at("x").get<double>() : 0.0,
                json.contains("y") ? json.at("y").get<double>() : 0.0,
                json.contains("z") ? json.at("z").get<double>() : 0.0});
    }


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
        IMatterCPtr const wind{
            json.contains("safierWind") ? parseSafierWind(json.at("safierWind")) : nullptr
        };

        MatterTranslationCPtr const translation{
            json.contains("translation") ? parseTranslation(json.at("translation")) : nullptr
        };

        return std::make_shared<FlaredDisk const>(
            json.at("rInner").get<double>(),
            json.at("rOuter").get<double>(),
            json.at("rho0").get<double>(),
            json.at("h0").get<double>(),
            json.at("r0").get<double>(),
            json.at("alpha").get<double>(),
            json.at("beta").get<double>(),
            wind,
            translation);
    }


    IMatterCPtr parseSphereEnvelope(const nlohmann::json& json)
    {
        MatterTranslationCPtr const translation{
            json.contains("translation") ? parseTranslation(json.at("translation")) : nullptr
        };

        return std::make_shared<SphereEnvelope const>(
            json.at("rInner").get<double>(),
            json.at("rOuter").get<double>(),
            json.at("rho0").get<double>(),
            json.at("r0").get<double>(),
            json.at("alpha").get<double>(),
            translation);
    }


    IMatterCPtr parseGeometry(const nlohmann::json& json)
    {
        DATA_ASSERT(
            json.size() == 1,
            "geometry section of the json must contain exactly one top level element");

        if (json.contains("flaredDisk"))
        {
            return parseFlaredDisk(json.at("flaredDisk"));
        } else if (json.contains("sphereEnvelope")) {
            return parseSphereEnvelope(json.at("sphereEnvelope"));
        } else if (json.contains("max")) {
            nlohmann::json const& jsonList = json.at("max");

            DATA_ASSERT(
                !jsonList.empty(),
                "max section of the geometry must be nonempty");

            std::vector<IMatterCPtr> matterArray;
            matterArray.reserve(jsonList.size());

            for (const auto& matterJson : jsonList)
            {
                matterArray.emplace_back(parseGeometry(matterJson));
            }

            return std::make_shared<MatterArray const>(
                std::move(matterArray),
                MatterArray::max);
        } else if (json.contains("sum")) {
            nlohmann::json const& jsonList = json.at("sum");

            DATA_ASSERT(
                !jsonList.empty(),
                "sum section of the geometry must be nonempty");

            std::vector<IMatterCPtr> matterArray;
            matterArray.reserve(jsonList.size());

            for (const auto& matterJson : jsonList)
            {
                matterArray.emplace_back(parseGeometry(matterJson));
            }

            return std::make_shared<MatterArray const>(
                std::move(matterArray),
                MatterArray::sum);
        }

        return nullptr;
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
            json.at("nx").get<uint32_t>(),
            json.at("ny").get<uint32_t>(),
            json.at("nz").get<uint32_t>(),
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

        DATA_ASSERT(
            !observers->empty(),
            "observers section of the json must contain at least one observer");
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
    IMatterCPtr geometry = parseGeometry(j.at("geometry"));
    grid_ = parseGrid(j.at("grid"), dustJson.at("kappa").get<double>(), geometry);
    sources_ = parseSources(j.at("stars"), sourceParameters);
    parseObservers(observers, j.at("observers"));
}
