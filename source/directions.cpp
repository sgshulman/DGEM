#include <cmath>
#include <cstddef>
#include <cassert>
#include "directions.hpp"

namespace
{
    double const PI = 3.1415926;

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

            inline uint32_t midpointsNumber() const
            { return midpointsNumber_; }

            inline PointWithNeighbors const *neighbor(int const i) const
            { return neighbors_[i]; }

            inline PointWithNeighbors *midpoint(int const i)
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
            uint32_t midpointsNumber_{ 0 };
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

                for(size_t cnt=0; cnt!=dots_[i]->midpointsNumber(); ++cnt )
                {
                    if ( dots_[i]->neighbor(cnt) == dots_[j] )
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
            explicit IcosahedronMesh(uint32_t directionsLevels)
                : directionsNumber_{ 0 }
                , readyPoints_{ 0 }
                , nodes_{ new PointWithNeighbors[nodesNumber(directionsLevels)]{} }
                , triangles_{ nullptr }
            {
                buildInitialMesh();

                for(uint32_t cnt=1; cnt<directionsLevels; ++cnt)
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

                for (uint64_t cnt=0; cnt<directionsNumber_; ++cnt)
                {
                    points[cnt] = triangles_[cnt].median();
                }

                return points;
            }

            double* weights() const
            {
                auto weights = new double[ directionsNumber_ ];

                for (uint64_t cnt=0; cnt<directionsNumber_; ++cnt)
                {
                    weights[cnt] = triangles_[cnt].square() * directionsNumber_ / 4. / PI;
                }

                return weights;
            }

            uint64_t directionsNumber() const
            {
                return directionsNumber_;
            }

        private:
            static uint32_t const ICOSAHEDRON_FACES{ 20 };
            static uint32_t const ICOSAHEDRON_VERTICES{ 12 };

            static uint64_t nodesNumber(uint32_t const directionsLevels)
            {
                uint64_t directionPointsNumber = ICOSAHEDRON_VERTICES;
                uint64_t directionsNumber = ICOSAHEDRON_FACES;

                for (size_t cnt=1; cnt<directionsLevels; ++cnt)
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
                    nodes_[cnt+1].set(Vector3d{ std::cos(angle*PI/180),  std::sin(angle*PI/180),  0.5 });
                    nodes_[cnt+6].set(Vector3d{ std::cos((36.+angle)*PI/180),  std::sin((36.+angle)*PI/180),  -0.5});
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

            PointWithNeighbors *getMidpoint(uint64_t const triangleId, int const i, int const j)
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

                for(size_t cnt = 0; cnt < directionsNumber_; ++cnt)
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

                for (uint64_t cnt=0; cnt!=readyPoints_; ++cnt)
                {
                    nodes_[cnt].removeMidpoints();
                }

                directionsNumber_ *= 4;
            }

            uint64_t directionsNumber_;
            uint64_t readyPoints_;
            PointWithNeighbors *nodes_;
            SphericalTriangle *triangles_;
    };
}


// directions
Directions::Directions(uint32_t NumOfDirectionsLevels)
{
    IcosahedronMesh mesh(NumOfDirectionsLevels);
    points_ = mesh.points();
    w_	  = mesh.weights();
    directionsNumber_ = mesh.directionsNumber();
}
