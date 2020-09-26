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
    char const sDust[] = "dust";
    char const sGeometry[] = "geometry";
    char const sGrid[] = "grid";
    char const sStars[] = "stars";
    char const sObservers[] = "observers";
    char const sFlaredDisk[] = "flaredDisk";
    char const sTranslation[] = "translation";
    char const sSphereEnvelope[] = "sphereEnvelope";
    char const sFractalCloud[] = "fractalCloud";

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


    double extract_double(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.is_number_float() && !json.is_number_integer())
        {
            throw std::invalid_argument(std::string("Item ") + name + " from section " + section + " should be float.");
        }

        return json.get<double>();
    }


    double get_double(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.contains(name))
        {
            throw std::invalid_argument(std::string("Section ") + section + " should contain float item " + name + ".");
        }

        return extract_double(json.at(name), section, name);
    }


    double get_optional_double(const nlohmann::json& json, char const* section, char const* name, double defaultValue)
    {
        return json.contains(name) ? extract_double(json.at(name), section, name) : defaultValue;
    }


    std::int32_t get_int32(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.contains(name))
        {
            throw std::invalid_argument(std::string("Section ") + section + " should contain integer item " + name + ".");
        }

        auto const& item = json.at(name);
        if (!item.is_number_integer())
        {
            throw std::invalid_argument(std::string("Item ") + name + " from section " + section + " should be integer.");
        }

        return item.get<std::int32_t>();
    }


    std::uint32_t get_uint32(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.contains(name))
        {
            throw std::invalid_argument(std::string("Section ") + section + " should contain unsigned integer item " + name + ".");
        }

        auto const& item = json.at(name);
        if (!item.is_number_unsigned() && !(item.is_number_integer() && item.get<int64_t>() >= 0))
        {
            throw std::invalid_argument(std::string("Item ") + name + " from section " + section + " should be unsigned integer.");
        }

        return item.get<std::uint32_t>();
    }


    std::uint64_t get_uint64(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.contains(name))
        {
            throw std::invalid_argument(std::string("Section ") + section + " should contain unsigned integer item " + name + ".");
        }

        auto const& item = json.at(name);
        if (!item.is_number_unsigned() && !(item.is_number_integer() && item.get<int64_t>() >= 0))
        {
            throw std::invalid_argument(std::string("Item ") + name + " from section " + section + " should be unsigned integer.");
        }

        return item.get<std::uint64_t>();
    }


    MatterTranslationCPtr parseTranslation(const nlohmann::json& json)
    {
        checkParameters(json, sTranslation, {"precession", "nutation", "intrinsicRotation", "x", "y", "z"});

        return std::make_shared<MatterTranslation const>(
            get_optional_double(json, sTranslation, "precession", 0.0),
            get_optional_double(json, sTranslation, "nutation", 0.0),
            get_optional_double(json, sTranslation, "intrinsicRotation", 0.0),
            Vector3d{
                get_optional_double(json, sTranslation, "x", 0.0),
                get_optional_double(json, sTranslation, "y", 0.0),
                get_optional_double(json, sTranslation, "z", 0.0)});
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
        checkParameters(
            json,
            sFlaredDisk,
            {"safierWind", sTranslation, "hump", "rInner", "rOuter", "rho0", "h0", "r0", "alpha", "beta"});

        IMatterCPtr const wind{
            json.contains("safierWind") ? parseSafierWind(json.at("safierWind")) : nullptr
        };

        MatterTranslationCPtr const translation{
            json.contains(sTranslation) ? parseTranslation(json.at(sTranslation)) : nullptr
        };

        IDiskHumpCPtr const hump{
            json.contains("hump") ? parseDiskHump(json.at("hump")) : nullptr
        };

        return std::make_shared<FlaredDisk const>(
            get_double(json, sFlaredDisk, "rInner"),
            get_double(json, sFlaredDisk, "rOuter"),
            get_double(json, sFlaredDisk, "rho0"),
            get_double(json, sFlaredDisk, "h0"),
            get_double(json, sFlaredDisk, "r0"),
            get_double(json, sFlaredDisk, "alpha"),
            get_double(json, sFlaredDisk, "beta"),
            wind,
            translation,
            hump);
    }


    IMatterCPtr parseSphereEnvelope(const nlohmann::json& json)
    {
        checkParameters(json, sSphereEnvelope, {sTranslation, "rInner", "rOuter", "rho0", "r0", "alpha"});

        MatterTranslationCPtr const translation{
            json.contains(sTranslation) ? parseTranslation(json.at(sTranslation)) : nullptr
        };

        return std::make_shared<SphereEnvelope const>(
            get_double(json, sSphereEnvelope, "rInner"),
            get_double(json, sSphereEnvelope, "rOuter"),
            get_double(json, sSphereEnvelope, "rho0"),
            get_double(json, sSphereEnvelope, "r0"),
            get_double(json, sSphereEnvelope, "alpha"),
            translation);
    }


    IMatterCPtr parseFractalCloud(const nlohmann::json& json)
    {
        checkParameters(json, sFractalCloud, {"n", "max", "dCube", "rho0", "dotsN", "seed"});

        return std::make_shared<FractalCloud const>(
            get_uint32(json, sFractalCloud, "n"),
            get_double(json, sFractalCloud, "max"),
            get_double(json, sFractalCloud, "dCube"),
            get_double(json, sFractalCloud, "rho0"),
            get_uint32(json, sFractalCloud, "dotsN"),
            get_int32(json, sFractalCloud, "seed"));
    }


    IMatterCPtr parseGeometry(const nlohmann::json& json)
    {
        DATA_ASSERT(
            json.size() == 1,
            "geometry section of the json must contain exactly one top level element");

        checkParameters(json, sGeometry, {sFlaredDisk, sSphereEnvelope, sFractalCloud, "max", "sum"});

        if (json.contains(sFlaredDisk))
        {
            return parseFlaredDisk(json.at(sFlaredDisk));
        } else if (json.contains(sSphereEnvelope)) {
            return parseSphereEnvelope(json.at(sSphereEnvelope));
        } else if (json.contains(sFractalCloud)) {
            return parseFractalCloud(json.at(sFractalCloud));
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
        checkParameters(json, sDust, {"kappa", "albedo", "hgg", "pl", "pc", "sc"});

        return std::make_shared<Dust>(
            get_double(json, sDust, "albedo"),
            get_double(json, sDust, "hgg"),
            get_double(json, sDust, "pl"),
            get_double(json, sDust, "pc"),
            get_double(json, sDust, "sc"));
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
            get_uint32(json, cartesian, "nx"),
            get_uint32(json, cartesian, "ny"),
            get_uint32(json, cartesian, "nz"),
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
        checkParameters(json, sGrid, {"cartesian", "tetrahedral"});

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
        checkParameters(json, sObservers, {"rimage", "manual", "parallel", "meridian"});

        auto const rimage = get_double(json, sObservers, "rimage");

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
            std::uint32_t const numberOfObservers = get_uint32(parallellJson, observersParallel, "numberOfObservers");
            const auto viewTheta = get_double(parallellJson, observersParallel, "theta");
            observers->reserve(observers->size() + numberOfObservers);

            for (std::uint32_t i=0; i!=numberOfObservers; ++i)
            {
                observers->emplace_back(2*PI/numberOfObservers*i, radians(viewTheta), rimage);
            }
        }

        if (json.contains("meridian"))
        {
            nlohmann::json const& medianJson = json.at("meridian");
            char const observersMeridian[] = "observers::meridian";
            checkParameters(medianJson, observersMeridian, {"numberOfObservers", "phi"});
            std::uint32_t const numberOfObservers = get_uint32(medianJson, observersMeridian, "numberOfObservers");
            const auto viewPhi = get_double(medianJson, observersMeridian, "phi");
            observers->reserve(observers->size() + numberOfObservers);

            for (std::uint32_t i=0; i!=numberOfObservers; ++i)
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
            checkParameters(star, sStars, {"x", "y", "z", "l"});

            Vector3d const position{
                get_double(star, sStars, "x"),
                get_double(star, sStars, "y"),
                get_double(star, sStars, "z")};

            pointSources.emplace_back(
                position,
                grid->cellId(position),
                get_double(star, sStars, "l"));
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
        {"method parameters", sDust, sGeometry, sGrid, sStars, sObservers});

    SourceParameters sourceParameters{};
    char const methodParameters[] = "method parameters";
    nlohmann::json const& methodJson = j.at(methodParameters);

    checkParameters(
        methodJson,
        methodParameters,
        {"fMonteCarlo", "nphotons", "PrimaryDirectionsLevel", "iseed", "taumin", "nscat",
         "SecondaryDirectionsLevel", "NumOfPrimaryScatterings", "NumOfSecondaryScatterings", "MonteCarloStart"});

    sourceParameters.useMonteCarlo_ = methodJson.at("fMonteCarlo").get<bool>();
    sourceParameters.num_photons_ = get_uint64(methodJson, methodParameters, "nphotons");
    sourceParameters.PrimaryDirectionsLevel_ = get_uint32(methodJson, methodParameters, "PrimaryDirectionsLevel");
    iseed_ = get_int32(methodJson, methodParameters, "iseed");

    fMonteCarlo_ = sourceParameters.useMonteCarlo_;
    taumin_ = get_double(methodJson, methodParameters, "taumin");
    nscat_ = get_uint32(methodJson, methodParameters, "nscat");
    SecondaryDirectionsLevel_ = get_uint32(methodJson, methodParameters, "SecondaryDirectionsLevel");
    NumOfPrimaryScatterings_ = get_uint32(methodJson, methodParameters, "NumOfPrimaryScatterings");
    NumOfSecondaryScatterings_ = get_uint32(methodJson, methodParameters, "NumOfSecondaryScatterings");
    MonteCarloStart_ = get_uint32(methodJson, methodParameters, "MonteCarloStart");

    nlohmann::json const& dustJson = j.at(sDust);
    dust_ = parseDust(dustJson);
    IMatterCPtr geometry = parseGeometry(j.at(sGeometry));
    grid_ = parseGrid(j.at(sGrid), get_double(dustJson, sDust, "kappa"), geometry);
    sources_ = parseSources(j.at(sStars), sourceParameters, grid_);
    parseObservers(observers, j.at(sObservers));
}
