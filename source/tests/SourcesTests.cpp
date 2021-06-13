#include "../third-party/catch2/catch.hpp"
#include "../Observer.hpp"
#include "../Sources.hpp"
#include "../IGrid.hpp"
#include "../Random.hpp"

namespace
{
    class TestGrid : public IGrid
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

        int movePhotonAtDepth(Photon& /*ph*/, double /*tau*/, double /*tauold*/) const override
        {
            return 0;
        }

        int movePhotonAtRandomDepth(Photon& /*ph*/, Random */*ran*/) const override
        {
            return 0;
        }

        void peeloff(Photon /*ph*/, Observer& /*observer*/, IDustCRef /*dust*/) const override
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

        bool inside(const Vector3d& /*position*/) const override
        {
            return true;
        }
    };
}


TEST_CASE("Sphere source. Random", "[sources]")
{
    std::vector<PointSource> pointSources;
    std::vector<SphereSource> sphereSources;
    sphereSources.emplace_back(Vector3d{0.,0.,0.}, 0, 1., 1.);

    SourceParameters const sourceParameters{true, false, 10, 1};
    Sources sources(sourceParameters, std::move(pointSources), std::move(sphereSources));

    IGridCPtr grid = std::make_shared<TestGrid>();
    Random rand(-1);

    for (int i=0; i!=10; ++i)
    {
        Photon ph = sources.emitRandomPhoton(grid, &rand);
        REQUIRE(ph.pos() * ph.dir().vector() >= 0.);
    }
}


TEST_CASE("Star disc", "[sources]")
{
    SourceParameters const sourceParameters{true, false, 10, 1};

    std::vector<SphereSource> sphereSources;
    sphereSources.emplace_back(Vector3d{0.,0.,0.}, 0, 1., 1.);
    Sources sources(sourceParameters, std::vector<PointSource>(), std::move(sphereSources));

    IGridCPtr grid = std::make_shared<TestGrid>();

    std::vector<Observer> observers;
    observers.emplace_back(0, 0, 1.);

    sources.directPhotons(grid, &observers);

    int const BELT_NUMBER{ 10 };
    double f[BELT_NUMBER]{};

    for (std::uint64_t x=0; x!=200; ++x)
    {
        for (std::uint64_t y=0; y != 200; ++y)
        {
            double dx = 0.01 * static_cast<double>(x - 100);
            double dy = 0.01 * static_cast<double>(y - 100);

            double const r = std::sqrt(dx * dx + dy * dy);

            if (0. < r && r < 1.)
            {
                int id = std::min(BELT_NUMBER - 1, static_cast<int>(BELT_NUMBER * r));
                f[id] += observers[0].totalLuminosity(x, y) / (2 * PI * r);
            }
        }
    }

    double const minF = *std::min_element(f, f + BELT_NUMBER);
    double const maxF = *std::max_element(f, f + BELT_NUMBER);

    REQUIRE((maxF / minF - 1) < 0.2);
}
