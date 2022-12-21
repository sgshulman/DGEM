#include <cmath>
#include <cstddef>
#include <cassert>
#include "DebugUtils.hpp"
#include "Directions.hpp"
#include "MathUtils.hpp"

namespace
{
    // Class presents mesh nodes on the union sphere.
    // Holds the dot and it neighbors for creating more detailed mesh
    class PointWithNeighbors
    {
        public:
            PointWithNeighbors()
                : pos_{}
            {};

            void set(Vector3d const& position)
            {
                pos_ = position.normalized();
            }

            inline Vector3d vector() const
            {
                return pos_;
            }

            void setMiddle(PointWithNeighbors const *d1, PointWithNeighbors const *d2)
            {
                pos_ = (d1->pos_ + d2->pos_).normalized();
                midpointsNumber_ = 0;
            }

            inline std::uint32_t midpointsNumber() const
            { return midpointsNumber_; }

            inline PointWithNeighbors const *neighbor(std::uint32_t const i) const
            { return neighbors_[i]; }

            inline PointWithNeighbors *midpoint(std::uint32_t const i)
            { return midpoints_[i]; }

            void addMidpoint(PointWithNeighbors const *neighbor, PointWithNeighbors *midpoint)
            {
                neighbors_[midpointsNumber_] = neighbor;
                midpoints_[midpointsNumber_] = midpoint;
                ++midpointsNumber_;
            }

            inline void removeMidpoints()
            { midpointsNumber_ = 0; }

        private:
            Vector3d pos_;
            std::uint32_t midpointsNumber_{ 0 };
            PointWithNeighbors const *neighbors_[6]{};
            PointWithNeighbors *midpoints_[6]{};
    };


    // Element of a mesh on the union sphere
    class SphericalTriangle
    {
        public:
            SphericalTriangle() = default;

            void set(PointWithNeighbors *dot1, PointWithNeighbors *dot2, PointWithNeighbors *dot3)
            {
                dots_[0] = dot1;
                dots_[1] = dot2;
                dots_[2] = dot3;
            }

            Vector3d median() const
            {
                return (dots_[0]->vector() + dots_[1]->vector() + dots_[2]->vector()).normalized();
            }

            double square() const
            {
                double const triple = tripleProduct(
                    dots_[0]->vector(),
                    dots_[1]->vector(),
                    dots_[2]->vector());

                double const r0r1 = dots_[0]->vector() * dots_[1]->vector();
                double const r1r2 = dots_[1]->vector() * dots_[2]->vector();
                double const r0r2 = dots_[0]->vector() * dots_[2]->vector();

                return 2. * std::atan(triple / (1. + r0r1 + r1r2 + r0r2));
            }

            PointWithNeighbors *readyMeadpoint(int const i, int const j)
            {
                assert(0 <= i && i <= 2 && 0 <= j && j <= 2);

                for(std::uint32_t cnt=0; cnt!=dots_[i]->midpointsNumber(); ++cnt )
                {
                    if (dots_[i]->neighbor(cnt) == dots_[j])
                    {
                        return dots_[i]->midpoint(cnt);
                    }
                }

                return nullptr;
            }

            PointWithNeighbors *operator[](int i)
            { return dots_[i]; }

        private:
            PointWithNeighbors *dots_[3]{};
    };


    class IcosahedronMesh
    {
        public:
            explicit IcosahedronMesh(std::uint32_t directionsLevels)
                : directionsNumber_{ 0 }
                , readyPoints_{ 0 }
                , nodes_{ new PointWithNeighbors[nodesNumber(directionsLevels)]{} }
                , triangles_{ nullptr }
            {
                buildInitialMesh();

                for(std::uint32_t cnt=1; cnt<directionsLevels; ++cnt)
                {
                    splitMesh();
                }
            }

            ~IcosahedronMesh()
            {
                delete[] nodes_;
                delete[] triangles_;
            }

            Vector3d* points() const
            {
                auto points = new Vector3d[ directionsNumber_ ];

                for (std::uint64_t cnt=0; cnt<directionsNumber_; ++cnt)
                {
                    points[cnt] = triangles_[cnt].median();
                }

                return points;
            }

            double* weights() const
            {
                auto weights = new double[ directionsNumber_ ];

                for (std::uint64_t cnt=0; cnt<directionsNumber_; ++cnt)
                {
                    weights[cnt] = triangles_[cnt].square() * static_cast<double>(directionsNumber_) / 4. / PI;
                }

                return weights;
            }

            std::uint64_t directionsNumber() const
            {
                return directionsNumber_;
            }

        private:
            static std::uint32_t const ICOSAHEDRON_FACES{ 20 };
            static std::uint32_t const ICOSAHEDRON_VERTICES{ 12 };

            static std::uint64_t nodesNumber(std::uint32_t const directionsLevels)
            {
                std::uint64_t directionPointsNumber = ICOSAHEDRON_VERTICES;
                std::uint64_t directionsNumber = ICOSAHEDRON_FACES;

                for (std::uint64_t cnt=1; cnt<directionsLevels; ++cnt)
                {
                    directionPointsNumber = directionPointsNumber + (directionsNumber * 3) / 2;
                    directionsNumber *= 4;
                }

                return directionPointsNumber;
            }

            void buildInitialMesh()
            {
                directionsNumber_ = ICOSAHEDRON_FACES;
                readyPoints_ = ICOSAHEDRON_VERTICES;
                triangles_ = new SphericalTriangle[directionsNumber_]{};

                // initial sphere conditions
                nodes_[0].set(Vector3d{0.0, 0.0, std::sqrt(5.0)/2.0});
                for (int cnt=0; cnt!=5; ++cnt)
                {
                    double const angle = cnt * 72.;
                    nodes_[cnt+1].set(Vector3d{ std::cos(radians(angle)),  std::sin(radians(angle)),  0.5 });
                    nodes_[cnt+6].set(Vector3d{ std::cos(radians(36.+angle)),  std::sin(radians(36.+angle)),  -0.5});
                }
                nodes_[11].set(Vector3d{0.0, 0.0, -std::sqrt(5.0)/2.0});

                // triangles
                triangles_[0].set(&nodes_[0], &nodes_[1], &nodes_[2]);
                triangles_[1].set(&nodes_[0], &nodes_[2], &nodes_[3]);
                triangles_[2].set(&nodes_[0], &nodes_[3], &nodes_[4]);
                triangles_[3].set(&nodes_[0], &nodes_[4], &nodes_[5]);
                triangles_[4].set(&nodes_[0], &nodes_[5], &nodes_[1]);
                triangles_[5].set(&nodes_[10], &nodes_[11], &nodes_[6]);
                triangles_[6].set(&nodes_[11], &nodes_[7], &nodes_[6]);
                triangles_[7].set(&nodes_[11], &nodes_[8], &nodes_[7]);
                triangles_[8].set(&nodes_[11], &nodes_[9], &nodes_[8]);
                triangles_[9].set(&nodes_[11], &nodes_[10], &nodes_[9]);
                triangles_[10].set(&nodes_[6], &nodes_[1], &nodes_[10]);
                triangles_[11].set(&nodes_[2], &nodes_[6], &nodes_[7]);
                triangles_[12].set(&nodes_[3], &nodes_[7], &nodes_[8]);
                triangles_[13].set(&nodes_[4], &nodes_[8], &nodes_[9]);
                triangles_[14].set(&nodes_[5], &nodes_[9], &nodes_[10]);
                triangles_[15].set(&nodes_[6], &nodes_[2], &nodes_[1]);
                triangles_[16].set(&nodes_[7], &nodes_[3], &nodes_[2]);
                triangles_[17].set(&nodes_[8], &nodes_[4], &nodes_[3]);
                triangles_[18].set(&nodes_[9], &nodes_[5], &nodes_[4]);
                triangles_[19].set(&nodes_[10], &nodes_[1], &nodes_[5]);
            }

            PointWithNeighbors *getMidpoint(std::uint64_t const triangleId, int const i, int const j)
            {
                assert(0 <= i && i <= 2 && 0 <= j && j <= 2);

                // search for ready midpoint
                PointWithNeighbors *midpoint = triangles_[triangleId].readyMeadpoint(i, j);

                if (!midpoint)
                {
                    nodes_[readyPoints_].setMiddle(triangles_[triangleId][i], triangles_[triangleId][j]);
                    midpoint = &(nodes_[readyPoints_]);
                    ++readyPoints_;
                    triangles_[triangleId][i]->addMidpoint(triangles_[triangleId][j], midpoint);
                    triangles_[triangleId][j]->addMidpoint(triangles_[triangleId][i], midpoint);
                }

                return midpoint;
            }

            void splitMesh()
            {
                auto newTriangles = new SphericalTriangle[directionsNumber_ * 4]{};

                for(std::uint64_t cnt = 0; cnt < directionsNumber_; ++cnt)
                {
                    PointWithNeighbors *d_rib01 = getMidpoint(cnt, 0, 1);
                    PointWithNeighbors *d_rib02 = getMidpoint(cnt, 0, 2);
                    PointWithNeighbors *d_rib12 = getMidpoint(cnt, 1, 2);

                    newTriangles[ 4 * cnt     ].set(triangles_[cnt][0], d_rib01, d_rib02);
                    newTriangles[ 4 * cnt + 1 ].set(triangles_[cnt][1], d_rib12, d_rib01);
                    newTriangles[ 4 * cnt + 2 ].set(triangles_[cnt][2], d_rib02, d_rib12);
                    newTriangles[ 4 * cnt + 3 ].set(d_rib01, d_rib12, d_rib02);
                }

                delete[] triangles_;
                triangles_ = newTriangles;

                for (std::uint64_t cnt=0; cnt!=readyPoints_; ++cnt)
                {
                    nodes_[cnt].removeMidpoints();
                }

                directionsNumber_ *= 4;
            }

            std::uint64_t directionsNumber_;
            std::uint64_t readyPoints_;
            PointWithNeighbors *nodes_;
            SphericalTriangle *triangles_;
    };

    std::pair<std::int64_t, std::int64_t> unpackNestedId(std::uint64_t squareId0)
    {
        std::uint64_t const squareId1 = (squareId0 & 0x9999999999999999U) | ((squareId0 & 0x4444444444444444U) >> 1U) | ((squareId0 & 0x2222222222222222U) << 1U);
        std::uint64_t const squareId2 = (squareId1 & 0xC3C3C3C3C3C3C3C3U) | ((squareId1 & 0x3030303030303030U) >> 2U) | ((squareId1 & 0x0C0C0C0C0C0C0C0CU) << 2U);
        std::uint64_t const squareId3 = (squareId2 & 0xF00FF00FF00FF00FU) | ((squareId2 & 0x0F000F000F000F00U) >> 4U) | ((squareId2 & 0x00F000F000F000F0U) << 4U);
        std::uint64_t const squareId4 = (squareId3 & 0xFF0000FFFF0000FFU) | ((squareId3 & 0x00FF000000FF0000U) >> 8U) | ((squareId3 & 0x0000FF000000FF00U) << 8U);
        std::uint64_t const squareId5 = (squareId4 & 0xFFFF00000000FFFFU) | ((squareId4 & 0x0000FFFF00000000U) >>16U) | ((squareId4 & 0x00000000FFFF0000U) <<16U);

        std::int64_t const x = static_cast<int64_t>(squareId5 & 0x00000000FFFFFFFFU);
        std::int64_t const y = static_cast<int64_t>((squareId5 & 0xFFFFFFFF00000000U) >> 32U);

        return {x, y};
    }
} // namespace


// directions
Directions::Directions(std::uint32_t const NumOfDirectionsLevels, bool const useHEALPixGrid)
{
    if (NumOfDirectionsLevels > 0)
    {
        if (useHEALPixGrid)
        {
            isolatitudeGrid3(4, NumOfDirectionsLevels);
        } else {
            IcosahedronMesh mesh(NumOfDirectionsLevels);
            points_ = mesh.points();
            w_ = mesh.weights();
            directionsNumber_ = mesh.directionsNumber();
        }
    }
}


// HEALPix like equal-area isolatitude tessellations of the sphere with N_theta = 3
Directions::Directions(std::uint32_t const Ntheta, std::uint32_t const Nphi, std::uint32_t const NumOfDirectionsLevels, bool const ringScheme)
{
    if (NumOfDirectionsLevels > 0)
    {
        DATA_ASSERT(Ntheta >= 2 && Ntheta <= 3, "Ntheta can be 2 or 3");
        DATA_ASSERT(Nphi % 2 == 0, "Nphi must be even");

        if (Ntheta == 3)
        {
            if (ringScheme)
            {
                isolatitudeGrid3(Nphi, NumOfDirectionsLevels);
            } else { // Nested Scheme
                isolatitudeGrid3Nested(Nphi, NumOfDirectionsLevels);
            }
        }
        else
        {
            if (ringScheme)
            {
                isolatitudeGrid2(Nphi, NumOfDirectionsLevels);
            } else { // Nested Scheme
                isolatitudeGrid2Nested(Nphi, NumOfDirectionsLevels);
            }
        }
    }
}


void Directions::isolatitudeGrid3(std::uint32_t Nphi, std::uint32_t const Nside)
{
    assert(Nphi % 2 == 0);

    directionsNumber_ = 3 * Nphi * Nside * Nside;
    points_ = new Vector3d[directionsNumber_]{};

    std::uint32_t p = 0;

    // Polar caps
    for (std::uint32_t i = 1, j = 1; p!= Nphi * Nside * (Nside - 1) / 2; ++j, ++p)
    {
        if (j > Nphi * i)
        {
            j = 1;
            ++i;
        }

        assert(1 <= i && i < Nside);
        assert(1 <= j && j <= i * Nphi);

        double const z = 1 - (i*i) / (3. * Nside * Nside);
        double const phi = 2 * PI / (Nphi * i) * (j - 0.5);
        points_[p] = Vector3d(phi, std::acos(z));
        points_[directionsNumber_ - 1 - p] = Vector3d(2 * PI - phi, std::acos(-z));
    }

    // Equatorial belts
    for (; p!= 3 * Nphi * Nside * Nside / 2; ++p)
    {
        std::uint32_t const pp = p - Nphi * Nside * (Nside - 1) / 2;
        std::uint32_t const i = pp / (Nphi * Nside) + Nside;
        std::uint32_t const j = pp % (Nphi * Nside) + 1;

        assert(Nside <= i && i <= 2 * Nside);
        assert(1 <= j && j <= Nphi * Nside);

        double const z = 4.0 / 3 - (2 * i) / (3. * Nside);
        std::uint32_t const s = (i - Nside) % 2;

        double const phi = 2 * PI / (Nphi * Nside) * (j - 0.5 * (1 + s));
        double const phi2 = 2 * PI / (Nphi * Nside) * (j - 0.5 * (1 - s));

        points_[p] = Vector3d(phi, std::acos(z));
        points_[directionsNumber_ - 1 - p] = Vector3d(2 * PI - phi2, std::acos(-z));
    }

    w_ = new double[directionsNumber_]{};

    for (std::uint32_t i=0; i!=directionsNumber_; ++i)
    {
        w_[i] = 1.0;
    }
}


void Directions::isolatitudeGrid3Nested(std::uint32_t Nphi, std::uint32_t const Nside)
{
    directionsNumber_ = 3 * Nphi * Nside * Nside;
    points_ = new Vector3d[directionsNumber_]{};

    for (std::uint64_t id = 0; id != directionsNumber_; ++id)
    {
        std::uint64_t const f = id / (Nside * Nside);
        std::uint64_t const frow = f / Nphi;

        std::uint64_t const squareId0 = id % (Nside * Nside);
        std::pair<std::int64_t, std::int64_t> xy = unpackNestedId(squareId0);
        std::int64_t const v = xy.first + xy.second;
        std::int64_t const h = xy.first - xy.second;

        std::int64_t const F1 = static_cast<std::int64_t>(frow) + 2;
        std::int64_t const F2 = static_cast<std::int64_t>(2 * (f % Nphi) - (frow % 2) + 1);

        std::int64_t i = F1 * Nside - v - 1;
        double z = 4.0 / 3.0 - 2.0 * static_cast<double>(i) / (3.0 * Nside);
        std::int64_t s = (i - Nside + 1) % 2;
        std::int64_t r = Nside;

        if (i < Nside) // North polar cap
        {
            r = i;
            z = 1.0 - static_cast<double>(i * i) / (3.0 * Nside * Nside);
            s = 1;
        }
        else if (i > 3 * Nside) // South polar cap
        {
            r = 4 * Nside - i;
            z = -1.0 + static_cast<double>(r * r) / (3.0 * Nside * Nside);
            s = 1;
        }

        std::int64_t j = (F2 * r + h + s) / 2;

        if (j > Nphi * Nside)
        {
            j -= Nphi * Nside;
        }

        if (j < 1)
        {
            j += Nphi * Nside;
        }

        points_[id] = Vector3d((static_cast<double>(j) - 0.5 * static_cast<double>(s)) * 2 * PI / static_cast<double>(Nphi * r), std::acos(z));
    }

     w_ = new double[directionsNumber_]{};

     for (std::uint32_t i=0; i!=directionsNumber_; ++i)
     {
        w_[i] = 1.0;
     }
}


void Directions::isolatitudeGrid2(std::uint32_t Nphi, std::uint32_t const Nside)
{
    directionsNumber_ = 2 * Nphi * Nside * Nside;
    points_ = new Vector3d[directionsNumber_]{};

    std::uint32_t p = 0;

    // Polar caps
    for (std::uint32_t i = 1, j = 1; p!= Nphi * Nside * (Nside - 1) / 2; ++j, ++p)
    {
        if (j > Nphi * i)
        {
            j = 1;
            ++i;
        }

        assert(1 <= i && i <= Nside);
        assert(1 <= j && j <= i * Nphi);

        double const z = 1 - (i*i) / (2. * Nside * Nside);
        double const phi = 2 * PI / (Nphi * i) * (j - 0.5);
        double const phi2 = 2 * PI / (Nphi * i) * j;
        points_[p] = Vector3d(phi, std::acos(z));
        points_[directionsNumber_ - 1 - p] = Vector3d(2 * PI - phi2, std::acos(-z));
    }

    // Equatorial belts
    for (; p!= Nphi * Nside * (3 * Nside + 1) / 2; ++p)
    {
        std::uint32_t const pp = p - Nphi * Nside * (Nside - 1) / 2;
        std::uint32_t const i = pp / (Nphi * Nside) + Nside;
        std::uint32_t const j = pp % (Nphi * Nside) + 1;

        assert(Nside <= i && i <= 2 * Nside);
        assert(1 <= j && j <= Nphi * Nside);

        double const z = 3.0 / 2 - (1. * i) / Nside;
        std::uint32_t const s = (i - Nside) % 2;

        double const phi = 2 * PI / (Nphi * Nside) * (j - 0.5 * (1 + s));

        points_[p] = Vector3d(phi, std::acos(z));
    }

    w_ = new double[directionsNumber_]{};

    for (std::uint32_t i=0; i!=directionsNumber_; ++i)
    {
        w_[i] = 1.0;
    }
}


void Directions::isolatitudeGrid2Nested(std::uint32_t Nphi, std::uint32_t const Nside)
{
    directionsNumber_ = 2 * Nphi * Nside * Nside;
    points_ = new Vector3d[directionsNumber_]{};

    for (std::uint64_t id = 0; id != directionsNumber_; ++id)
    {
        std::uint64_t const f = id / (Nside * Nside);
        std::uint64_t const frow = f / Nphi;

        std::uint64_t const squareId0 = id % (Nside * Nside);
        std::pair<std::int64_t, std::int64_t> xy = unpackNestedId(squareId0);
        std::int64_t const v = xy.first + xy.second;
        std::int64_t const h = xy.first - xy.second;

        std::int64_t const F1 = static_cast<std::int64_t>(frow) + 2;
        std::int64_t const F2 = static_cast<std::int64_t>(2 * (f % Nphi) - (frow % 2) + 1);

        std::int64_t i = F1 * Nside - v - 1;
        double z = 3.0 / 2 - static_cast<double>(i) / Nside;
        std::int64_t s = (i - Nside + 1) % 2;
        std::int64_t r = Nside;

        if (i < Nside) // North polar cap
        {
            r = i;
            z = 1 - static_cast<double>(r * r) / (2. * Nside * Nside);
            s = 1;
        }
        else if (i > 2 * Nside) // South polar cap
        {
            r = 3 * Nside - i;
            z = -1.0 + static_cast<double>(r * r) / (2.0 * Nside * Nside);
            s = 0;
        }

        std::int64_t j = (F2 * r + h + s) / 2;

        if (j > Nphi * Nside)
        {
            j -= Nphi * Nside;
        }

        if (j < 1)
        {
            j += Nphi * Nside;
        }

        points_[id] = Vector3d((static_cast<double>(j) - 0.5 * static_cast<double>(s)) * 2 * PI / static_cast<double>(Nphi * r), std::acos(z));
    }

     w_ = new double[directionsNumber_]{};

     for (std::uint32_t i=0; i!=directionsNumber_; ++i)
     {
        w_[i] = 1.0;
     }
}
