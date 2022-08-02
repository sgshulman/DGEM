#include <iostream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include "model.hpp"
#include "AzimuthalHump.hpp"
#include "CartesianGrid.hpp"
#include "DebugUtils.hpp"
#include "WhiteDust.hpp"
#include "FlaredDisk.hpp"
#include "FractalCloud.hpp"
#include "KurosawaWind.hpp"
#include "MathUtils.hpp"
#include "MatterArray.hpp"
#include "MatterTranslation.hpp"
#include "MieDust.hpp"
#include "LEcuyer.hpp"
#include "RoundHump.hpp"
#include "SafierWind.hpp"
#include "Sobol.hpp"
#include "Sources.hpp"
#include "SphereEnvelope.hpp"
#include "TetrahedralGrid.hpp"
#include "Observer.hpp"
#include "third-party/nlohmann/json.hpp"
#include "Units.hpp"

namespace
{
    char const sDensitySlice[] = "densitySlice";
    char const sDust[] = "dust";
    char const sEffectiveHeight[] = "effectiveHeight";
    char const sGeometry[] = "geometry";
    char const sGrid[] = "grid";
    char const sStars[] = "stars";
    char const sObservers[] = "observers";
    char const sFlaredDisk[] = "flaredDisk";
    char const sSafierWind[] = "safierWind";
    char const sKurosawaWind[] = "kurosawaWind";
    char const sMeanOpticalDepth[] = "meanOpticalDepth";
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
                ss << "Section \"" << section << "\" contains invalid item " << k.key()
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


    std::string extract_string(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.is_string())
        {
            throw std::invalid_argument(std::string("Item ") + name + " from section " + section + " should be string.");
        }

        return json.get<std::string>();
    }


    std::string get_string(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.contains(name))
        {
            throw std::invalid_argument(std::string("Section ") + section + " should contain string item " + name + ".");
        }

        return extract_string(json.at(name), section, name);
    }


    std::string get_optional_string(const nlohmann::json& json, char const* section, char const* name, std::string const& defaultValue)
    {
        return json.contains(name) ? extract_string(json.at(name), section, name) : defaultValue;
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


    bool extract_bool(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.is_boolean())
        {
            throw std::invalid_argument(std::string("Item ") + name + " from section " + section + " should be boolean.");
        }

        return json.get<bool>();
    }

    bool get_bool(const nlohmann::json& json, char const* const section, char const* const name)
    {
        if (!json.contains(name))
        {
            throw std::invalid_argument(std::string("Section ") + section + " should contain boolean item " + name + ".");
        }

        return extract_bool(json.at(name), section, name);
    }

    bool get_optional_bool(const nlohmann::json& json, char const* section, char const* name, bool defaultValue)
    {
        return json.contains(name) ? extract_bool(json.at(name), section, name) : defaultValue;
    }

    template<typename T>
    T extract_enum(
        const nlohmann::json& json,
        char const* section,
        char const* name,
        std::initializer_list<std::string> enumValues)
    {
        if (!json.is_string())
        {
            throw std::invalid_argument(std::string("Item ") + name + " from section " + section + " should be string.");
        }

        std::string const str = json.get<std::string>();
        std::uint32_t i = 0;
        for (auto const& value : enumValues)
        {
            if (str == value)
            {
                return static_cast<T>(i);
            }
            ++i;
        }

        std::stringstream ss;
        ss << "Enum \"" << name << "\" from section \"" << section << "\" contains invalid value " << str
           << ". Possible values are: ";

        for (const auto& value : enumValues)
        {
            ss << value << ", ";
        }
        ss << "\n";

        throw std::invalid_argument(ss.str());

        return T();
    }

    template<typename T>
    T get_optional_enum(
        const nlohmann::json& json,
        char const* section,
        char const* name,
        std::initializer_list<std::string> enumValues,
        T defaultValue)
    {
        return json.contains(name) ? extract_enum<T>(json.at(name), section, name, enumValues) : defaultValue;
    }


    MatterTranslationCPtr parseTranslation(const nlohmann::json& json)
    {
        checkParameters(json, sTranslation, {"precession", "nutation", "intrinsicRotation", "x", "y", "z"});

        return std::make_shared<MatterTranslation const>(
            radians(get_optional_double(json, sTranslation, "precession", 0.0)),
            radians(get_optional_double(json, sTranslation, "nutation", 0.0)),
            radians(get_optional_double(json, sTranslation, "intrinsicRotation", 0.0)),
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
            checkParameters(jsonHump, azimuthalHump, {"h", "r", "sigma2", "sigma2azimuthal", "sigma2azimuthalBackward"});

            return std::make_shared<AzimuthalHump const>(
                get_double(jsonHump, azimuthalHump, "h"),
                get_double(jsonHump, azimuthalHump, "r"),
                get_double(jsonHump, azimuthalHump, "sigma2"),
                get_double(jsonHump, azimuthalHump, "sigma2azimuthal"),
                get_optional_double(jsonHump, azimuthalHump, "sigma2azimuthalBackward", 0.0));
        }

        return nullptr;
    }


    IMatterCPtr parseSafierWind(const nlohmann::json& json)
    {
        checkParameters(json, sSafierWind, {"hump", "model", "mOut", "mStar", "h0", "rMin", "rMax"});

        IDiskHumpCPtr const hump{
            json.contains("hump") ? parseDiskHump(json.at("hump")) : nullptr
        };

        std::string const windModel = get_string(json, sSafierWind, "model");
        DATA_ASSERT(windModel.size() == 1, "Safier wind model should contain one letter.");

        return std::make_shared<SafierWind const>(
            windModel.at(0),
            get_double(json, sSafierWind, "mOut"),
            get_double(json, sSafierWind, "mStar"),
            get_double(json, sSafierWind, "h0"),
            get_double(json, sSafierWind, "rMin"),
            get_double(json, sSafierWind, "rMax"),
            hump);
    }


    IMatterCPtr parseKurosawaWind(const nlohmann::json& json)
    {
        checkParameters(
            json,
            sKurosawaWind,
            {"d", "p", "mOut", "mStar", "rInner", "rOuter", "rScale", "windAccelerationRate", "terminalV", "soundSpeed"});

        return std::make_shared<KurosawaWind const>(
            get_double(json, sKurosawaWind, "d"),
            get_double(json, sKurosawaWind, "p"),
            get_double(json, sKurosawaWind, "mOut"),
            get_double(json, sKurosawaWind, "mStar"),
            get_double(json, sKurosawaWind, "rInner"),
            get_double(json, sKurosawaWind, "rOuter"),
            get_double(json, sKurosawaWind, "rScale"),
            get_double(json, sKurosawaWind, "windAccelerationRate"),
            get_double(json, sKurosawaWind, "terminalV"),
            get_double(json, sKurosawaWind, "soundSpeed"));
    }


    IMatterCPtr parseFlaredDisk(const nlohmann::json& json)
    {
        checkParameters(
            json,
            sFlaredDisk,
            {sSafierWind, sKurosawaWind, sTranslation, "hump", "rInner", "rOuter", "rho0", "h0", "r0", "alpha", "beta"});

        IMatterCPtr const safierWind = json.contains(sSafierWind) ? parseSafierWind(json.at(sSafierWind)) : nullptr;
        IMatterCPtr const kurosawaWind = json.contains(sKurosawaWind) ? parseKurosawaWind(json.at(sKurosawaWind)) : nullptr;

        IMatterCPtr wind;

        if (safierWind && kurosawaWind)
        {
            wind = std::make_shared<MatterArray const>(
                std::vector<IMatterCPtr>{safierWind, kurosawaWind},
                MatterArray::max);

        } else if (safierWind) {
            wind = safierWind;
        } else if (kurosawaWind) {
            wind = kurosawaWind;
        }

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


    IDustCPtr parseWhiteDust(const nlohmann::json& json)
    {
        char const white[] = "dust::white";
        checkParameters(json, white, {"albedo", "hgg", "pl", "pc", "sc"});

        return std::make_shared<WhiteDust>(
            get_double(json, white, "albedo"),
            get_double(json, white, "hgg"),
            get_double(json, white, "pl"),
            get_double(json, white, "pc"),
            get_double(json, white, "sc"));
    }


    IDustCPtr parseMieDust(const nlohmann::json& json)
    {
        char const mie[] = "dust::mie";
        checkParameters(json, mie, {"albedo", "tableFile"});

        return std::make_shared<MieDust>(
            get_double(json, mie, "albedo"),
            get_string(json, mie, "tableFile"));
    }


    IDustCPtr parseDust(const nlohmann::json& json)
    {
        checkParameters(json, sDust, {"kappa", "white", "mie"});

        DATA_ASSERT(
            json.size() == 2,
            "grid section of the json must contain exactly two top level elements");

        if (json.contains("white"))
        {
            return parseWhiteDust(json.at("white"));
        } else if (json.contains("mie")) {
            return parseMieDust(json.at("mie"));
        }

        return nullptr;
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

        DATA_ASSERT(
            json.contains("nodesFile") == json.contains("elementsFile"),
            "Grid::tetrahedral should contain both files nodesFile and elementsFile or none of them");

        if (json.contains("nodesFile") && json.contains("elementsFile"))
        {
            return std::make_shared<TetrahedralGrid const>(
                get_string(json, tetrahedral, "nodesFile"),
                get_string(json, tetrahedral, "elementsFile"),
                get_optional_string(json, tetrahedral, "gridBinFile", ""),
                get_double(json, tetrahedral, "max"),
                kappa,
                std::move(matter));
        }

        return std::make_shared<TetrahedralGrid const>(
            get_string(json, tetrahedral, "gridBinFile"),
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

    void densitySlice(nlohmann::json const& json, IMatterCPtr const& matter)
    {
        checkParameters(json, sDensitySlice, {"filename", "phi", "radiusMax", "heightMax"});

        std::string const filename = get_string(json, sDensitySlice, "filename");
        double const phi = radians(get_double(json, sDensitySlice, "phi"));
        double const radiusMax = get_double(json, sDensitySlice, "radiusMax");
        double const heightMax = get_double(json, sDensitySlice, "heightMax");
        std::int32_t radiusN{ 400 };
        std::int32_t heightN{ 300 };

        double const dRadius = radiusMax / (radiusN - 1);
        double const dHeight = heightMax / (heightN - 1);
        double const cosPhi = std::cos(phi);
        double const sinPhi = std::sin(phi);

        std::ofstream file(filename);

        for (std::int32_t i = 0; i != heightN; ++i)
        {
            double const z = i * dHeight;

            for (std::int32_t j = 0; j != radiusN; ++j)
            {
                double const r = j * dRadius;
                double const x = r * cosPhi;
                double const y = r * sinPhi;

                file << matter->density({x, y, z}) << "\t";
            }
            file << std::endl;
        }
    }

    void effectiveHeight(nlohmann::json const& json, IMatterCPtr const& matter, double const kappa)
    {
        checkParameters(json, sEffectiveHeight, {"filename", "radiusMax", "heightMax", "dHeight", "dRadius", "printCoordinates"});

        std::string const filename = get_string(json, sEffectiveHeight, "filename");
        double const radiusMax = get_double(json, sEffectiveHeight, "radiusMax");
        double const heightMax = get_double(json, sEffectiveHeight, "heightMax");
        double const dHeight = get_double(json, sEffectiveHeight, "dHeight");
        double const dRadius = get_double(json, sEffectiveHeight, "dRadius");
        bool const printCoordinates = get_optional_bool(json, sEffectiveHeight, "printCoordinates", false);

        std::int32_t radiusN{ static_cast<int32_t>(2 * radiusMax / dRadius) + 1};
        std::int32_t heightN{ static_cast<int32_t>(heightMax / dHeight) + 1};

        std::ofstream file(filename);

        for (std::int32_t i = 0; i != radiusN; ++i)
        {
            double const x = -radiusMax + i * dRadius;

            for (std::int32_t j = 0; j != radiusN; ++j)
            {
                double const y = -radiusMax + j * dRadius;

                double z{ heightMax };
                double tau{ 0. };
                for (std::int32_t k = 0; k != heightN; ++k)
                {
                    tau += dHeight * matter->density({x, y, z}) * kappa * AU_Cm;

                    if (tau > 1.0)
                    {
                        break;
                    }

                    z -= dHeight;
                }

                if (printCoordinates)
                {
                    file << x << "\t" << y << "\t" << z << "\n";
                } else {
                    file << z << "\t";
                }
            }
            file << std::endl;
        }
    }

    void meanOpticalDepth(nlohmann::json const& json, IGridCRef grid)
    {
        checkParameters(json, sMeanOpticalDepth, {"useHEALPixGrid", "directionsLevel"});

        bool const useHEALPixGrid = get_optional_bool(json, sMeanOpticalDepth, "useHEALPixGrid", false);
        std::uint32_t const directionsLevel = get_uint32(json, sMeanOpticalDepth, "directionsLevel");

        Directions directions(directionsLevel, useHEALPixGrid);
        double opticalDepth{ 0.0 };

        for (std::uint64_t i=0; i!=directions.number(); ++i)
        {
            opticalDepth += grid->findRealOpticalDepth(Vector3d{0, 0, 0}, directions.direction(i));
        }

        opticalDepth /= directions.number();
        std::cout << "MeanOpticalDepth: " << opticalDepth << std::endl;
    }
}

#ifdef ENABLE_UNIT_TESTS
#include "tests/ModelTests.inl"
#endif

Model::Model(std::vector<Observer>* observers, std::string const& parametersFile)
{
    std::ifstream configurationFile(parametersFile);
    DATA_ASSERT(configurationFile.is_open(), parametersFile + " should exist.");
    nlohmann::json j;
    configurationFile >> j;

    char const methodParameters[] = "method parameters";

    checkParameters(
        j,
        "parameters.json",
        {methodParameters, sDensitySlice, sDust, sEffectiveHeight, sGeometry, sGrid, sStars, sObservers, sMeanOpticalDepth});

    SourceParameters sourceParameters{};
    nlohmann::json const& methodJson = j.at(methodParameters);

    checkParameters(
        methodJson,
        methodParameters,
        {"fMonteCarlo", "nphotons", "PrimaryDirectionsLevel", "generatorType", "iseed", "dgemBinType", "inputRandomFile",
         "outputRandomFile", "taumin", "nscat", "SecondaryDirectionsLevel", "NumOfPrimaryScatterings",
         "NumOfSecondaryScatterings", "MonteCarloStart", "fUseHEALPixGrid", "defaultStarRadius", "fWriteScatterings"});

    sourceParameters.useMonteCarlo_ = get_bool(methodJson, methodParameters, "fMonteCarlo");
    sourceParameters.num_photons_ = get_uint64(methodJson, methodParameters, "nphotons");
    sourceParameters.PrimaryDirectionsLevel_ = get_uint32(methodJson, methodParameters, "PrimaryDirectionsLevel");
    sourceParameters.useHEALPixGrid_ = get_optional_bool(methodJson, methodParameters, "fUseHEALPixGrid", false);

    generatorType_ = get_optional_enum(methodJson, methodParameters, "generatorType", {"LEcuyer", "Sobol"}, RandomGeneratorType::LECUYER);
    if (LECUYER == generatorType_)
    {
        iseed_ = get_int32(methodJson, methodParameters, "iseed");
    }

    dgemBinType_ = get_optional_enum(methodJson, methodParameters, "dgemBinType", {"Point", "Line", "HexLines"}, DgemBinType::LINE);
    inputRandomFile_ = get_optional_string(methodJson, methodParameters, "inputRandomFile", "");
    outputRandomFile_ = get_optional_string(methodJson, methodParameters, "outputRandomFile", "");

    fMonteCarlo_ = sourceParameters.useMonteCarlo_;
    writeScatterings_ = get_optional_bool(methodJson, methodParameters, "fWriteScatterings", true);
    taumin_ = get_double(methodJson, methodParameters, "taumin");
    defaultStarRadius_ = get_optional_double(methodJson, methodParameters, "defaultStarRadius", 0.0047);
    nscat_ = get_uint32(methodJson, methodParameters, "nscat");
    SecondaryDirectionsLevel_ = get_uint32(methodJson, methodParameters, "SecondaryDirectionsLevel");
    NumOfPrimaryScatterings_ = get_uint32(methodJson, methodParameters, "NumOfPrimaryScatterings");
    NumOfSecondaryScatterings_ = get_uint32(methodJson, methodParameters, "NumOfSecondaryScatterings");
    MonteCarloStart_ = get_uint32(methodJson, methodParameters, "MonteCarloStart");
    useHEALPixGrid_ = sourceParameters.useHEALPixGrid_;

    nlohmann::json const& dustJson = j.at(sDust);
    dust_ = parseDust(dustJson);
    IMatterCPtr geometry = parseGeometry(j.at(sGeometry));

    if (j.contains(sDensitySlice))
    {
        densitySlice(j.at(sDensitySlice), geometry);
    }

    grid_ = parseGrid(j.at(sGrid), get_double(dustJson, sDust, "kappa"), geometry);

    if (j.contains(sEffectiveHeight))
    {
        effectiveHeight(j.at(sEffectiveHeight), geometry, get_double(dustJson, sDust, "kappa"));
    }

    if (j.contains(sMeanOpticalDepth))
    {
        meanOpticalDepth(j.at(sMeanOpticalDepth), grid_);
    }

    sources_ = parseSources(j.at(sStars), sourceParameters, grid_);
    parseObservers(observers, j.at(sObservers));
}

IRandomGenerator* Model::createRandomGenerator() const
{
    IRandomGenerator* rand{ nullptr };
    if (RandomGeneratorType::LECUYER == generatorType_)
    {
        rand = new LEcuyer(iseed_);
    }
    else if (RandomGeneratorType::SOBOL == generatorType_)
    {
        rand = new Sobol(3);
    }

    assert(rand != nullptr);

    if (!inputRandomFile_.empty())
    {
        rand->load(inputRandomFile_);
    }

    rand->setOutputFile(outputRandomFile_);
    return rand;
}
