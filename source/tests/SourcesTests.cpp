#include "../third-party/catch2/catch.hpp"
#include "../Sources.hpp"
#include "../IGrid.hpp"
#include "../LEcuyer.hpp"

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

        bool movePhotonAtDepth(Photon& /*ph*/, double /*tau*/, double /*tauold*/) const override
        {
            return 0;
        }

        bool movePhotonAtRandomDepth(Photon& /*ph*/, IRandomGenerator */*ran*/) const override
        {
            return 0;
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
}


TEST_CASE("Sphere source. Random", "[sources]")
{
    std::vector<PointSource> pointSources;
    std::vector<SphereSource> sphereSources;
    sphereSources.emplace_back(Vector3d{0.,0.,0.}, 0, 1., 1.);

    SourceParameters const sourceParameters{10, 1, true, false};
    Sources sources(sourceParameters, std::move(pointSources), std::move(sphereSources));

    IGridCPtr grid = std::make_shared<TestGrid>();
    LEcuyer rand(-1);

    for (int i=0; i!=10; ++i)
    {
        Photon ph = sources.emitRandomPhoton(grid, &rand);
        REQUIRE(ph.pos() * ph.dir().vector() >= 0.);
    }
}
