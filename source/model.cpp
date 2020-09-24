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

    double get_double(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.contains(name))
        {
            std::stringstream ss;
            ss << "Section " << section << " should contain a float type item " << name << ".";

            throw std::invalid_argument(ss.str());
        }
        const auto item = json.at(name);

        if (!item.is_number_float() && !item.is_number_integer())
        {
            std::stringstream ss;
            ss << "Item " << name << " from section " << section << " should be float.";

            throw std::invalid_argument(ss.str());
        }

        return item.get<double>();
    }


    double get_optional_double(
        const nlohmann::json& json,
        char const* const section,
        char const* const name,
        double defaultValue)
    {
        if (!json.contains(name))
        {
            return defaultValue;
        }

        const auto item = json.at(name);

        if (!item.is_number_float() && !item.is_number_integer())
        {
            std::stringstream ss;
            ss << "Item " << name << " from section " << section << " should be float.";

            throw std::invalid_argument(ss.str());
        }

        return item.get<double>();
    }


    MatterTranslationCPtr parseTranslation(const nlohmann::json& json)
    {
        const char transition[] = "translation";
        checkParameters(json, transition, {"precession", "nutation", "intrinsicRotation", "x", "y", "z"});

        return std::make_shared<MatterTranslation const>(
            get_optional_double(json, transition, "precession", 0.0),
            get_optional_double(json, transition, "nutation", 0.0),
            get_optional_double(json, transition, "intrinsicRotation", 0.0),
            Vector3d{
                get_optional_double(json, transition, "x", 0.0),
                get_optional_double(json, transition, "y", 0.0),
                get_optional_double(json, transition, "z", 0.0)});
    }

    IDiskHumpCPtr parseDiskHump(const nlohmann::json& json)
    {
        if (json.empty())
        {
            return nullptr;
        }

        DATA_ASSERT(json.size() == 1, "only one hump for the element is allowed");

        char const roundHump[] = "roundHump";
        char const azimuthalHump[] = "azimuthalHump";
        checkParameters(json, "hump", {roundHump, azimuthalHump});

        if (json.contains("roundHump"))
        {
            nlohmann::json const& jsonHump = json.at(roundHump);
            checkParameters(jsonHump, roundHump, {"h", "r", "sigma2"});

            return std::make_shared<RoundHump const>(
                get_double(jsonHump, roundHump, "h"),
                get_double(jsonHump, roundHump, "r"),
                get_double(jsonHump, roundHump, "sigma2"));
        } else if (json.contains("azimuthalHump")) {
            nlohmann::json const& jsonHump = json.at(azimuthalHump);
            checkParameters(jsonHump, azimuthalHump, {"h", "r", "sigma2", "sigma2azimuthal"});

            return std::make_shared<AzimuthalHump const>(
                get_double(jsonHump, azimuthalHump, "h"),
                get_double(jsonHump, azimuthalHump, "r"),
                get_double(jsonHump, azimuthalHump, "sigma2"),
                get_double(jsonHump, azimuthalHump, "sigma2azimuthal"));
        }

        return nullptr;
    }


    IMatterCPtr parseSafierWind(const nlohmann::json& json)
    {
        char const safierWind[] = "safierWind";
        checkParameters(json, safierWind, {"hump", "model", "mOut", "mStar", "h0", "rMin", "rMax"});

        IDiskHumpCPtr const hump{
            json.contains("hump") ? parseDiskHump(json.at("hump")) : nullptr
        };

        return std::make_shared<SafierWind const>(
            json.at("model").get<std::string>().at(0),
            get_double(json, safierWind, "mOut"),
            get_double(json, safierWind, "mStar"),
            get_double(json, safierWind, "h0"),
            get_double(json, safierWind, "rMin"),
            get_double(json, safierWind, "rMax"),
            hump);
    }


    IMatterCPtr parseFlaredDisk(const nlohmann::json& json)
    {
        char const flaredDisk[] = "flaredDisk";

        checkParameters(
            json,
            flaredDisk,
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
            get_double(json, flaredDisk, "rInner"),
            get_double(json, flaredDisk, "rOuter"),
            get_double(json, flaredDisk, "rho0"),
            get_double(json, flaredDisk, "h0"),
            get_double(json, flaredDisk, "r0"),
            get_double(json, flaredDisk, "alpha"),
            get_double(json, flaredDisk, "beta"),
            wind,
            translation,
            hump);
    }


    IMatterCPtr parseSphereEnvelope(const nlohmann::json& json)
    {
        char const sphereEnvelope[] = "sphereEnvelope";
        checkParameters(json, sphereEnvelope, {"translation", "rInner", "rOuter", "rho0", "r0", "alpha"});

        MatterTranslationCPtr const translation{
            json.contains("translation") ? parseTranslation(json.at("translation")) : nullptr
        };

        return std::make_shared<SphereEnvelope const>(
            get_double(json, sphereEnvelope, "rInner"),
            get_double(json, sphereEnvelope, "rOuter"),
            get_double(json, sphereEnvelope, "rho0"),
            get_double(json, sphereEnvelope, "r0"),
            get_double(json, sphereEnvelope, "alpha"),
            translation);
    }


    IMatterCPtr parseFractalCloud(const nlohmann::json& json)
    {
        char const fractalCloud[] = "fractalCloud";
        checkParameters(json, fractalCloud, {"n", "max", "dCube", "rho0", "dotsN", "seed"});

        return std::make_shared<FractalCloud const>(
            json.at("n").get<std::uint32_t>(),
            get_double(json, fractalCloud, "max"),
            get_double(json, fractalCloud, "dCube"),
            get_double(json, fractalCloud, "rho0"),
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
        char const dust[] = "dust";
        checkParameters(json, dust, {"kappa", "albedo", "hgg", "pl", "pc", "sc"});

        return std::make_shared<Dust>(
            get_double(json, dust, "albedo"),
            get_double(json, dust, "hgg"),
            get_double(json, dust, "pl"),
            get_double(json, dust, "pc"),
            get_double(json, dust, "sc"));
    }


    IGridCPtr parseCartesianGrid(const nlohmann::json& json, double const kappa, IMatterCPtr matter)
    {
        char const cartesian[] = "grid::cartesian";
        checkParameters(json, cartesian, {"xmax", "ymax", "zmax", "nx", "ny", "nz"});

        return std::make_shared<CartesianGrid const>(
            get_double(json, cartesian, "xmax"),
            get_double(json, cartesian, "ymax"),
            get_double(json, cartesian, "zmax"),
            kappa,
            json.at("nx").get<std::uint32_t>(),
            json.at("ny").get<std::uint32_t>(),
            json.at("nz").get<std::uint32_t>(),
            std::move(matter));
    }


    IGridCPtr parseTetrahedralGrid(const nlohmann::json& json, double const kappa, IMatterCPtr matter)
    {
        char const tetrahedral[] = "grid::tetrahedral";
        checkParameters(json, "tetrahedral", {"nodesFile", "elementsFile", "gridBinFile", "max"});

        if (json.contains("nodesFile") && json.contains("elementsFile"))
        {
            return std::make_shared<TetrahedralGrid const>(
                json.at("nodesFile").get<std::string>(),
                json.at("elementsFile").get<std::string>(),
                json.contains("gridBinFile") ? json.at("gridBinFile").get<std::string>() : "",
                get_double(json, tetrahedral, "max"),
                kappa,
                std::move(matter));
        }

        return std::make_shared<TetrahedralGrid const>(
            json.at("gridBinFile").get<std::string>(),
            get_double(json, tetrahedral, "max"),
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
        char const name[] = "observers";
        checkParameters(json, name, {"rimage", "manual", "parallel", "meridian"});

        auto const rimage = get_double(json, name, "rimage");

        if (json.contains("manual"))
        {
            nlohmann::json const& manualJson = json.at("manual");
            observers->reserve(manualJson.size());
            char const observersManual[] = "observers::manual";
            for (const auto& observer : manualJson)
            {
                checkParameters(observer, observersManual, {"phi", "theta"});

                observers->emplace_back(
                    radians(get_double(observer, observersManual, "phi")),
                    radians(get_double(observer, observersManual, "theta")),
                    rimage);
            }
        }

        if (json.contains("parallel"))
        {
            nlohmann::json const& parallellJson = json.at("parallel");
            char const observersParallel[] = "observers::parallel";
            checkParameters(parallellJson, observersParallel, {"numberOfObservers", "theta"});
            const auto numberOfObservers = parallellJson.at("numberOfObservers").get<int>();
            const auto viewTheta = get_double(parallellJson, observersParallel, "theta");
            observers->reserve(observers->size() + numberOfObservers);

            for (int i=0; i!=numberOfObservers; ++i)
            {
                observers->emplace_back(2*PI/numberOfObservers*i, radians(viewTheta), rimage);
            }
        }

        if (json.contains("meridian"))
        {
            nlohmann::json const& medianJson = json.at("meridian");
            char const observersMeridian[] = "observers::meridian";
            checkParameters(medianJson, observersMeridian, {"numberOfObservers", "phi"});
            const auto numberOfObservers = medianJson.at("numberOfObservers").get<int>();
            const auto viewPhi = get_double(medianJson, observersMeridian, "phi");
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
        char const sources[] = "sources";

        for (auto const &star : json)
        {
            checkParameters(star, sources, {"x", "y", "z", "l"});

            Vector3d const position{
                get_double(star, sources, "x"),
                get_double(star, sources, "y"),
                get_double(star, sources, "z")};

            pointSources.emplace_back(
                position,
                grid->cellId(position),
                get_double(star, sources, "l"));
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
    char const methodParameters[] = "method parameters";
    nlohmann::json const& methodJson = j.at(methodParameters);

    checkParameters(
        methodJson,
        methodParameters,
        {"fMonteCarlo", "nphotons", "PrimaryDirectionsLevel", "iseed", "taumin", "nscat",
         "SecondaryDirectionsLevel", "NumOfPrimaryScatterings", "NumOfSecondaryScatterings", "MonteCarloStart"});

    sourceParameters.useMonteCarlo_ = methodJson.at("fMonteCarlo").get<bool>();
    sourceParameters.num_photons_ = methodJson.at("nphotons").get<std::uint64_t>();
    sourceParameters.PrimaryDirectionsLevel_ = methodJson.at("PrimaryDirectionsLevel").get<std::uint32_t>();
    iseed_ = methodJson.at("iseed").get<int32_t>();

    fMonteCarlo_ = sourceParameters.useMonteCarlo_;
    taumin_ = get_double(methodJson, methodParameters, "taumin");
    nscat_ = methodJson.at("nscat").get<std::uint32_t>();
    SecondaryDirectionsLevel_ = methodJson.at("SecondaryDirectionsLevel").get<std::uint32_t>();
    NumOfPrimaryScatterings_ = methodJson.at("NumOfPrimaryScatterings").get<std::uint32_t>();
    NumOfSecondaryScatterings_ = methodJson.at("NumOfSecondaryScatterings").get<std::uint32_t>();
    MonteCarloStart_ = methodJson.at("MonteCarloStart").get<std::uint32_t>();

    char const dust[] = "dust";
    nlohmann::json const& dustJson = j.at(dust);
    dust_ = parseDust(dustJson);
    IMatterCPtr geometry = parseGeometry(j.at("geometry"));
    grid_ = parseGrid(j.at("grid"), get_double(dustJson, dust, "kappa"), geometry);
    sources_ = parseSources(j.at("stars"), sourceParameters, grid_);
    parseObservers(observers, j.at("observers"));
}
