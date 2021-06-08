#include "../third-party/catch2/catch.hpp"
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
