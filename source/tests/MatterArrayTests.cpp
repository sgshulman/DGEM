#include "../third-party/catch2/catch.hpp"
#include "../MatterArray.hpp"

namespace
{
    class MatterStub : public IMatter
    {
        public:
            explicit MatterStub(double rho)
                : rho_{ rho }
            {}

            ~MatterStub() override = default;

            double density(Vector3d const& /*position*/) const override
            {
                return rho_;
            }

        private:
            double const rho_;
    };
} // namespace

TEST_CASE("Matter Array", "[matter]")
{
    SECTION("Maximum")
    {
        MatterArray array{
            std::vector<IMatterCPtr>{
                std::make_shared<MatterStub>(1.),
                std::make_shared<MatterStub>(2.),
                std::make_shared<MatterStub>(3.)},
            MatterArray::max};

        REQUIRE(Approx(3.) == array.density({0., 0., 0.}));
    }

    SECTION("Sum")
    {
        MatterArray array{
            std::vector<IMatterCPtr>{
                std::make_shared<MatterStub>(1.),
                std::make_shared<MatterStub>(2.),
                std::make_shared<MatterStub>(3.)},
            MatterArray::sum};

        REQUIRE(Approx(6.) == array.density({0., 0., 0.}));
    }
}
