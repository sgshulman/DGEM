#include "../third-party/catch2/catch.hpp"
#include "../Observer.hpp"
#include "../Sources.hpp"
#include "../IGrid.hpp"
#include "../LEcuyer.hpp"

#include <sstream>

namespace
{
    class TestGrid final: public IGrid
    {
        double findRealOpticalDepth(Vector3d const& /*position*/, Vector3d const& /*direction*/) const override
        {
            return 0.;
        }

        double findOpticalDepth(Photon /*ph*/) const override
        {
            return 0.;
        }

        double movePhotonAtDistance(Photon& ph, double distance) const override
        {
            ph.Move(distance, 0);
            return 0.;
        }

        bool movePhotonAtDepth(Photon& /*ph*/, double /*tau*/, double /*tauold*/) const override
        {
            return true;
        }

        bool movePhotonAtRandomDepth(Photon& /*ph*/, IRandomGenerator* /*ran*/) const override
        {
            return true;
        }

        void peeloff(Photon /*ph*/, Observer& /*observer*/, IDustCRef /*dust*/) const override
        {}

        void peeloff(Photon /*ph*/, Observer& /*observer*/, IDustCRef /*dust*/, const Vector3d& /*pos1*/, const Vector3d& /*pos2*/) const override
        {}

        void peeloffHex(Photon /*ph*/, Observer& /*observer*/, IDustCRef /*dust*/, Vector3d const& /*pos1*/, Vector3d const& /*pos2*/) const override
        {}

        double computeMatterMass() const override
        {
            return 0.;
        }

        double max() const override
        {
            return 100;
        }

        std::uint64_t cellId(const Vector3d& /*position*/) const override
        {
            return 0;
        }

        bool inside(const Photon& /*ph*/) const override
        {
            return true;
        }
    };

    double parseTotalLuminosity(std::stringstream& log)
    {
        std::string phiStr;
        double phi;
        std::string thetaStr;
        double theta;
        std::string fStr;
        double f;

        log >> phiStr >> phi >> thetaStr >> theta >> fStr >> f;
        return f;
    }

    double sourceLuminosity(std::uint64_t num_photons, double radius, std::uint32_t diskLevel)
    {
        SourceParameters const sourceParameters{num_photons, 0, 0, diskLevel, true, false};

        std::vector<SphereSource> sphereSources;
        std::vector<PointSource> pointSources;

        if (radius < std::numeric_limits<double>::epsilon())
        {
            sphereSources.emplace_back(Vector3d{0., 0., 0.}, 0, 1., radius);
        } else {
            pointSources.emplace_back(Vector3d{0., 0., 0.}, 0, 1.);
        }

        Sources sources(sourceParameters, std::move(pointSources), std::move(sphereSources));
        IGridCPtr grid = std::make_shared<TestGrid>();

        std::vector<Observer> observers;
        observers.emplace_back(0, 0, 100., 0., 200, 200, 1);

        sources.directPhotons(grid, &observers);

        std::stringstream log;
        observers[0].normalize(num_photons);
        observers[0].write(log);

        return parseTotalLuminosity(log);
    }
} // namespace


TEST_CASE("Sphere source. Random", "[sources]")
{
    std::vector<PointSource> pointSources;
    std::vector<SphereSource> sphereSources;
    sphereSources.emplace_back(Vector3d{0.,0.,0.}, 0, 1., 1.);

    SourceParameters const sourceParameters{10, 0, 0, 100, true, false};
    Sources sources(sourceParameters, std::move(pointSources), std::move(sphereSources));

    IGridCPtr grid = std::make_shared<TestGrid>();
    LEcuyer rand(-1);

    for (int i=0; i!=10; ++i)
    {
        double randomValue;
        Photon ph = sources.emitRandomPhoton(grid, &rand, &randomValue);
        REQUIRE(ph.pos() * ph.dir().vector() >= 0.);
    }
}


TEST_CASE("Source Luminosity", "[sources]")
{
    double const sphereSourceLuminosity10 = sourceLuminosity(10, 1., 100);
    double const sphereSourceLuminosity100 = sourceLuminosity(100, 1., 100);
    double const sphereSourceLuminositySmall = sourceLuminosity(100, 0.1, 100);
    double const sphereSourceLuminosityLowLevel = sourceLuminosity(100, 1., 10);
    double const pointSourceLuminosity10 = sourceLuminosity(10, 0., 1);
    double const pointSourceLuminosity100 = sourceLuminosity(100, 0., 1);

    REQUIRE(Approx(pointSourceLuminosity10) == sphereSourceLuminosity10);
    REQUIRE(Approx(pointSourceLuminosity10) == pointSourceLuminosity100);
    REQUIRE(Approx(sphereSourceLuminosity10) == sphereSourceLuminosity100);
    REQUIRE(Approx(sphereSourceLuminosity10) == sphereSourceLuminositySmall);
    REQUIRE(Approx(sphereSourceLuminosity10) == sphereSourceLuminosityLowLevel);
}


TEST_CASE("Star disc", "[sources]")
{
    SourceParameters const sourceParameters{10, 0, 0, 128, true, false};

    std::vector<SphereSource> sphereSources;
    sphereSources.emplace_back(Vector3d{0.,0.,0.}, 0, 1., 1.);
    Sources sources(sourceParameters, std::vector<PointSource>(), std::move(sphereSources));

    IGridCPtr grid = std::make_shared<TestGrid>();

    std::vector<Observer> observers;
    observers.emplace_back(0, 0, 1., 0., 200, 200, 1);

    sources.directPhotons(grid, &observers);

    int const BELT_NUMBER{ 10 };
    double f[BELT_NUMBER]{};

    for (std::int32_t x=0; x!=200; ++x)
    {
        for (std::int32_t y=0; y != 200; ++y)
        {
            double dx = 0.01 * (1.0 * x - 99.5);
            double dy = 0.01 * (1.0 * y - 99.5);

            double const r = std::sqrt(dx * dx + dy * dy);

            int id = std::min(BELT_NUMBER - 1, static_cast<int>(BELT_NUMBER * r));
            f[id] += observers[0].totalLuminosity(x, y) / (2 * PI * r);
        }
    }

    double const minF = *std::min_element(f, f + BELT_NUMBER);
    double const maxF = *std::max_element(f, f + BELT_NUMBER);

    REQUIRE((maxF / minF - 1) < 0.05);
}

TEST_CASE("Sphere star is opaque for direct light", "[sources]")
{
    std::uint64_t const numPhotons = 100;
    SourceParameters const sourceParameters{numPhotons, 0, 0, 10, true, false};

    std::vector<SphereSource> sphereSources;
    std::vector<PointSource> pointSources;

    sphereSources.emplace_back(Vector3d{10., 0., 0.}, 0, 1., 1.);
    pointSources.emplace_back(Vector3d{-10., 0., 0.}, 0, 1.);

    Sources sources(sourceParameters, std::move(pointSources), std::move(sphereSources));
    IGridCPtr grid = std::make_shared<TestGrid>();

    std::vector<Observer> observers;
    observers.emplace_back(0, 0, 100., 0., 200, 200, 1);
    observers.emplace_back(0, radians(90), 100., 0., 200, 200, 1);
    observers.emplace_back(radians(90), radians(90), 100., 0., 200, 200, 1);
    observers.emplace_back(radians(180), radians(90), 100., 0., 200, 200, 1);

    sources.directPhotons(grid, &observers);

    std::stringstream log0;
    observers[0].normalize(numPhotons);
    observers[0].write(log0);

    double const poleLuminosity = parseTotalLuminosity(log0);

    std::stringstream log1;
    observers[1].normalize(numPhotons);
    observers[1].write(log1);
    double const equatorLuminosity = parseTotalLuminosity(log1);

    std::stringstream log2;
    observers[2].normalize(numPhotons);
    observers[2].write(log2);
    double const equatorLuminosity2 = parseTotalLuminosity(log2);

    std::stringstream log3;
    observers[3].normalize(numPhotons);
    observers[3].write(log3);
    double const equatorLuminosity3 = parseTotalLuminosity(log3);

    REQUIRE(Approx(poleLuminosity) == 1 / 4. / PI);
    REQUIRE(Approx(poleLuminosity) == 2. * equatorLuminosity);
    REQUIRE(Approx(poleLuminosity) == equatorLuminosity2);
    REQUIRE(Approx(poleLuminosity) == equatorLuminosity3);
}

TEST_CASE("intersectSphereSource", "[sources]")
{
    SourceParameters const sourceParameters{100, 0, 0, 10, true, false};

    std::vector<SphereSource> sphereSources;
    std::vector<PointSource> pointSources;
    sphereSources.emplace_back(Vector3d{0., 0., 0.}, 0, 1., 1.);

    Sources sources(sourceParameters, std::move(pointSources), std::move(sphereSources));

    REQUIRE(sources.intersectSphereSource({-1.05, 0., 0.}, {1., 0., 0}));
    REQUIRE(sources.intersectSphereSource({0., 0., 0.}, {1., 0., 0}));
    REQUIRE(!sources.intersectSphereSource({1.05, 0., 0.}, {1., 0., 0}));

    REQUIRE(sources.intersectSphereSource({0., 0.95, 0.}, {1., 0., 0}));
    REQUIRE(sources.intersectSphereSource({0., -0.95, 0.}, {1., 0., 0}));
    REQUIRE(sources.intersectSphereSource({0., 0., 0.95}, {1., 0., 0}));
    REQUIRE(sources.intersectSphereSource({0., 0., -0.95}, {1., 0., 0}));

    REQUIRE(!sources.intersectSphereSource({0., 1.05, 0.}, {1., 0., 0}));
    REQUIRE(!sources.intersectSphereSource({0., -1.05, 0.}, {1., 0., 0}));
    REQUIRE(!sources.intersectSphereSource({0., 0., 1.05}, {1., 0., 0}));
    REQUIRE(!sources.intersectSphereSource({0., 0., -1.05}, {1., 0., 0}));

    REQUIRE(!sources.intersectSphereSource({0.71, 0.71, 0.}, {1., 0., 0}));
    REQUIRE(!sources.intersectSphereSource({0.71, 0., -0.71}, {1., 0., 0}));
    REQUIRE(sources.intersectSphereSource({0.7, 0.7, 0.}, {1., 0., 0}));
    REQUIRE(sources.intersectSphereSource({0.7, 0., -0.7}, {1., 0., 0}));

    REQUIRE(sources.intersectSphereSource({0., -1.05, 0.}, {0., 1., 0}));
    REQUIRE(sources.intersectSphereSource({0., 0., 0.}, {0., 0., 0}));
    REQUIRE(!sources.intersectSphereSource({0., 1.05, 0.}, {0., 1., 0}));

    REQUIRE(sources.intersectSphereSource(Vector3d{1., 1., 1}.normalized() * 0.95,  Vector3d{1., 1., 1}.normalized()));
    REQUIRE(!sources.intersectSphereSource(Vector3d{1., 1., 1}.normalized() * 1.05,  Vector3d{1., 1., 1}.normalized()));

    REQUIRE(sources.intersectSphereSource(Vector3d{1., -1., 1}.normalized() * 0.95,  Vector3d{1., -1., 1}.normalized()));
    REQUIRE(!sources.intersectSphereSource(Vector3d{1., -1., 1}.normalized() * 1.05,  Vector3d{1., -1., 1}.normalized()));

    REQUIRE(sources.intersectSphereSource(Vector3d{-1., 1., -1}.normalized() * 0.95,  Vector3d{-1., 1., -1}.normalized()));
    REQUIRE(!sources.intersectSphereSource(Vector3d{-1., 1., -1}.normalized() * 1.05,  Vector3d{-1., 1., -1}.normalized()));
}

