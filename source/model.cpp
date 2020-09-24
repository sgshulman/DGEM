#include <iostream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include "model.hpp"
#include "AzimuthalHump.hpp"
#include "CartesianGrid.hpp"
#include "DebugUtils.hpp"
#include "Dust.hpp"
#include "FlaredDisk.hpp"
#include "FractalCloud.hpp"
#include "MathUtils.hpp"
#include "MatterArray.hpp"
#include "MatterTranslation.hpp"
#include "RoundHump.hpp"
#include "SafierWind.hpp"
#include "Sources.hpp"
#include "SphereEnvelope.hpp"
#include "TetrahedralGrid.hpp"
#include "observers.hpp"
#include "third-party/nlohmann/json.hpp"

namespace
{
    void checkParameters(
        const nlohmann::json& json,
        char const* const section,
        std::initializer_list<const char*> items)
    {
        for (auto k = json.cbegin(); k != json.cend(); ++k)
        {
            if (cend(items) == std::find(cbegin(items), cend(items), k.key()))
            {
                std::stringstream ss;
                ss << "Section " << section << " contains invalid item " << k.key()
                   << ". Possible items are: ";

                for (const auto& i: items)
                {
                    ss << i << ", ";
                }
                ss << "\n";

                throw std::invalid_argument(ss.str());
            }
        }
    }

    MatterTranslationCPtr parseTranslation(const nlohmann::json& json)
    {
        checkParameters(json, "translation", {"precession", "nutation", "intrinsicRotation", "x", "y", "z"});

        return std::make_shared<MatterTranslation const>(
            json.contains("precession") ? radians(json.at("precession").get<double>()) : 0.0,
            json.contains("nutation") ? radians(json.at("nutation").get<double>()) : 0.0,
            json.contains("intrinsicRotation") ? radians(json.at("intrinsicRotation").get<double>()) : 0.0,
            Vector3d{
                json.contains("x") ? json.at("x").get<double>() : 0.0,
                json.contains("y") ? json.at("y").get<double>() : 0.0,
                json.contains("z") ? json.at("z").get<double>() : 0.0});
    }

    IDiskHumpCPtr parseDiskHump(const nlohmann::json& json)
    {
        if (json.empty())
        {
            return nullptr;
        }

        DATA_ASSERT(json.size() == 1, "only one hump for the element is allowed");

        checkParameters(json, "hump", {"roundHump", "azimuthalHump"});

        if (json.contains("roundHump"))
        {
            nlohmann::json const& jsonHump = json.at("roundHump");
            checkParameters(jsonHump, "roundHump", {"h", "r", "sigma2"});

            return std::make_shared<RoundHump const>(
                jsonHump.at("h").get<double>(),
                jsonHump.at("r").get<double>(),
                jsonHump.at("sigma2").get<double>());
        } else if (json.contains("azimuthalHump")) {
            nlohmann::json const& jsonHump = json.at("azimuthalHump");
            checkParameters(jsonHump, "azimuthalHump", {"h", "r", "sigma2", "sigma2azimuthal"});

            return std::make_shared<AzimuthalHump const>(
                jsonHump.at("h").get<double>(),
                jsonHump.at("r").get<double>(),
                jsonHump.at("sigma2").get<double>(),
                jsonHump.at("sigma2azimuthal").get<double>());
        }

        return nullptr;
    }


    IMatterCPtr parseSafierWind(const nlohmann::json& json)
    {
        checkParameters(json, "safierWind", {"hump", "model", "mOut", "mStar", "h0", "rMin", "rMax"});

        IDiskHumpCPtr const hump{
            json.contains("hump") ? parseDiskHump(json.at("hump")) : nullptr
        };

        return std::make_shared<SafierWind const>(
            json.at("model").get<std::string>().at(0),
            json.at("mOut").get<double>(),
            json.at("mStar").get<double>(),
            json.at("h0").get<double>(),
            json.at("rMin").get<double>(),
            json.at("rMax").get<double>(),
            hump);
    }


    IMatterCPtr parseFlaredDisk(const nlohmann::json& json)
    {
        checkParameters(
            json,
            "flaredDisk",
            {"safierWind", "translation", "hump", "rInner", "rOuter", "rho0", "h0", "r0", "alpha", "beta"});

        IMatterCPtr const wind{
            json.contains("safierWind") ? parseSafierWind(json.at("safierWind")) : nullptr
        };

        MatterTranslationCPtr const translation{
            json.contains("translation") ? parseTranslation(json.at("translation")) : nullptr
        };

        IDiskHumpCPtr const hump{
            json.contains("hump") ? parseDiskHump(json.at("hump")) : nullptr
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
            translation,
            hump);
    }


    IMatterCPtr parseSphereEnvelope(const nlohmann::json& json)
    {
        checkParameters(json, "sphereEnvelope", {"translation", "rInner", "rOuter", "rho0", "r0", "alpha"});

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


    IMatterCPtr parseFractalCloud(const nlohmann::json& json)
    {
        checkParameters(json, "fractalCloud", {"n", "max", "dCube", "rho0", "dotsN", "seed"});

        return std::make_shared<FractalCloud const>(
            json.at("n").get<std::uint32_t>(),
            json.at("max").get<double>(),
            json.at("dCube").get<double>(),
            json.at("rho0").get<double>(),
            json.at("dotsN").get<std::uint32_t>(),
            json.at("seed").get<int32_t>());
    }


    IMatterCPtr parseGeometry(const nlohmann::json& json)
    {
        DATA_ASSERT(
            json.size() == 1,
            "geometry section of the json must contain exactly one top level element");

        checkParameters(json, "geometry", {"flaredDisk", "sphereEnvelope", "fractalCloud", "max", "sum"});

        if (json.contains("flaredDisk"))
        {
            return parseFlaredDisk(json.at("flaredDisk"));
        } else if (json.contains("sphereEnvelope")) {
            return parseSphereEnvelope(json.at("sphereEnvelope"));
        } else if (json.contains("fractalCloud")) {
            return parseFractalCloud(json.at("fractalCloud"));
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
        checkParameters(json, "dust", {"kappa", "albedo", "hgg", "pl", "pc", "sc"});

        return std::make_shared<Dust>(
            json.at("albedo").get<double>(),
            json.at("hgg").get<double>(),
            json.at("pl").get<double>(),
            json.at("pc").get<double>(),
            json.at("sc").get<double>());
    }


    IGridCPtr parseCartesianGrid(const nlohmann::json& json, double const kappa, IMatterCPtr matter)
    {
        checkParameters(json, "cartesian", {"xmax", "ymax", "zmax", "nx", "ny", "nz"});

        return std::make_shared<CartesianGrid const>(
            json.at("xmax").get<double>(),
            json.at("ymax").get<double>(),
            json.at("zmax").get<double>(),
            kappa,
            json.at("nx").get<std::uint32_t>(),
            json.at("ny").get<std::uint32_t>(),
            json.at("nz").get<std::uint32_t>(),
            std::move(matter));
    }


    IGridCPtr parseTetrahedralGrid(const nlohmann::json& json, double const kappa, IMatterCPtr matter)
    {
        checkParameters(json, "tetrahedral", {"nodesFile", "elementsFile", "gridBinFile", "max"});

        if (json.contains("nodesFile") && json.contains("elementsFile"))
        {
            return std::make_shared<TetrahedralGrid const>(
                json.at("nodesFile").get<std::string>(),
                json.at("elementsFile").get<std::string>(),
                json.contains("gridBinFile") ? json.at("gridBinFile").get<std::string>() : "",
                json.at("max").get<double>(),
                kappa,
                std::move(matter));
        }

        return std::make_shared<TetrahedralGrid const>(
            json.at("gridBinFile").get<std::string>(),
            json.at("max").get<double>(),
            kappa,
            std::move(matter));
    }


    IGridCPtr parseGrid(const nlohmann::json& json, double const kappa, IMatterCPtr matter)
    {
        checkParameters(json, "grid", {"cartesian", "tetrahedral"});

        DATA_ASSERT(
            json.size() == 1,
            "grid section of the json must contain exactly one top level element");

        if (json.contains("cartesian"))
        {
            return parseCartesianGrid(json.at("cartesian"), kappa, std::move(matter));
        } else if (json.contains("tetrahedral")) {
            return parseTetrahedralGrid(json.at("tetrahedral"), kappa, std::move(matter));
        }

        return nullptr;
    }


    void parseObservers(std::vector<Observer>* observers, nlohmann::json const& json)
    {
        checkParameters(json, "observers", {"rimage", "manual", "parallel", "meridian"});

        auto const rimage = json.at("rimage").get<double>();

        if (json.contains("manual"))
        {
            nlohmann::json const& manualJson = json.at("manual");
            observers->reserve(manualJson.size());
            for (const auto& observer : manualJson)
            {
                checkParameters(observer, "observers::manual", {"phi", "theta"});

                observers->emplace_back(
                    radians(observer.at("phi").get<double>()),
                    radians(observer.at("theta").get<double>()),
                    rimage);
            }
        }

        if (json.contains("parallel"))
        {
            nlohmann::json const& parallellJson = json.at("parallel");
            checkParameters(parallellJson, "observers::parallel", {"numberOfObservers", "theta"});
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
            checkParameters(medianJson, "observers::meridian", {"numberOfObservers", "phi"});
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

    SourcesPtr parseSources(nlohmann::json const& json, SourceParameters const& sourceParameters, IGridCRef grid)
    {
        std::vector<PointSource> pointSources;
        pointSources.reserve(json.size());

        for (auto const &star : json)
        {
            checkParameters(star, "sources", {"x", "y", "z", "l"});

            Vector3d const position{
                star.at("x").get<double>(),
                star.at("y").get<double>(),
                star.at("z").get<double>()};

            pointSources.emplace_back(
                position,
                grid->cellId(position),
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

    checkParameters(
        j,
        "parameters.json",
        {"method parameters", "dust", "geometry", "grid", "stars", "observers"});

    SourceParameters sourceParameters{};
    nlohmann::json const& methodJson = j.at("method parameters");

    checkParameters(
        methodJson,
        "method parameters",
        {"fMonteCarlo", "nphotons", "PrimaryDirectionsLevel", "iseed", "taumin", "nscat",
         "SecondaryDirectionsLevel", "NumOfPrimaryScatterings", "NumOfSecondaryScatterings", "MonteCarloStart"});

    sourceParameters.useMonteCarlo_ = methodJson.at("fMonteCarlo").get<bool>();
    sourceParameters.num_photons_ = methodJson.at("nphotons").get<std::uint64_t>();
    sourceParameters.PrimaryDirectionsLevel_ = methodJson.at("PrimaryDirectionsLevel").get<std::uint32_t>();
    iseed_ = methodJson.at("iseed").get<int32_t>();

    fMonteCarlo_ = sourceParameters.useMonteCarlo_;
    taumin_ = methodJson.at("taumin").get<double>();
    nscat_ = methodJson.at("nscat").get<std::uint32_t>();
    SecondaryDirectionsLevel_ = methodJson.at("SecondaryDirectionsLevel").get<std::uint32_t>();
    NumOfPrimaryScatterings_ = methodJson.at("NumOfPrimaryScatterings").get<std::uint32_t>();
    NumOfSecondaryScatterings_ = methodJson.at("NumOfSecondaryScatterings").get<std::uint32_t>();
    MonteCarloStart_ = methodJson.at("MonteCarloStart").get<std::uint32_t>();

    nlohmann::json const& dustJson = j.at("dust");
    dust_ = parseDust(dustJson);
    IMatterCPtr geometry = parseGeometry(j.at("geometry"));
    grid_ = parseGrid(j.at("grid"), dustJson.at("kappa").get<double>(), geometry);
    sources_ = parseSources(j.at("stars"), sourceParameters, grid_);
    parseObservers(observers, j.at("observers"));
}
