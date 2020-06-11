#include <cmath>
#include <cstddef>
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

            PointWithNeighbors *operator[](int i)
            { return dots_[i]; }

        private:
            PointWithNeighbors *dots_[3]{};
    };
}

// directions
Directions::Directions(uint32_t NumOfDirectionsLevels)
{
    NumOfDirections_ = 20;
    uint64_t NumOfReadyDirDots;
    uint64_t NumOfDirDots = 12;
    for(size_t cnt=1; cnt<NumOfDirectionsLevels; cnt++)
    {
        NumOfDirDots = NumOfDirDots+(NumOfDirections_*3)/2;
        NumOfDirections_ *= 4;
    }
    NumOfDirections_ = 20;

    auto *DirDots = new PointWithNeighbors[ NumOfDirDots ];
    auto *DirTrigon = new SphericalTriangle[ NumOfDirections_ ];
    SphericalTriangle *DirTrigon2;

    PointWithNeighbors *d_rib01, *d_rib02, *d_rib12;
    // initial sphere conditions
    (DirDots[0]).set( Vector3d{0.0, 0.0, std::sqrt(5.0)/2.0});
    for (int cnt=0; cnt!=5; ++cnt)
    {
        double const angle = cnt * 72.;
        DirDots[cnt+1].set(Vector3d{ std::cos(angle*PI/180),  std::sin(angle*PI/180),  0.5 });
        DirDots[cnt+6].set(Vector3d{ std::cos((36+angle)*PI/180),  std::sin((36+angle)*PI/180),  -0.5});
    }
    (DirDots[11]).set(Vector3d{0.0, 0.0, -std::sqrt(5.0)/2.0});
    NumOfReadyDirDots = 12;

    // triangles
    DirTrigon[0].set(&DirDots[0], &DirDots[1], &DirDots[2]);
    DirTrigon[1].set(&DirDots[0], &DirDots[2], &DirDots[3]);
    DirTrigon[2].set(&DirDots[0], &DirDots[3], &DirDots[4]);
    DirTrigon[3].set(&DirDots[0], &DirDots[4], &DirDots[5]);
    DirTrigon[4].set(&DirDots[0], &DirDots[5], &DirDots[1]);
    DirTrigon[5].set(&DirDots[10], &DirDots[11], &DirDots[6]);
    DirTrigon[6].set(&DirDots[11], &DirDots[7], &DirDots[6]);
    DirTrigon[7].set(&DirDots[11], &DirDots[8], &DirDots[7]);
    DirTrigon[8].set(&DirDots[11], &DirDots[9], &DirDots[8]);
    DirTrigon[9].set(&DirDots[11], &DirDots[10], &DirDots[9]);
    DirTrigon[10].set(&DirDots[6], &DirDots[1], &DirDots[10]);
    DirTrigon[11].set(&DirDots[2], &DirDots[6], &DirDots[7]);
    DirTrigon[12].set(&DirDots[3], &DirDots[7], &DirDots[8]);
    DirTrigon[13].set(&DirDots[4], &DirDots[8], &DirDots[9]);
    DirTrigon[14].set(&DirDots[5], &DirDots[9], &DirDots[10]);
    DirTrigon[15].set(&DirDots[6], &DirDots[2], &DirDots[1]);
    DirTrigon[16].set(&DirDots[7], &DirDots[3], &DirDots[2]);
    DirTrigon[17].set(&DirDots[8], &DirDots[4], &DirDots[3]);
    DirTrigon[18].set(&DirDots[9], &DirDots[5], &DirDots[4]);
    DirTrigon[19].set(&DirDots[10], &DirDots[1], &DirDots[5]);

    for(uint32_t cnt=1; cnt<NumOfDirectionsLevels; ++cnt)
    {
        // new triangles
        DirTrigon2 = DirTrigon;
        DirTrigon = new SphericalTriangle[NumOfDirections_ * 4];
        // splitting
        for(size_t cnt2 = 0; cnt2 < NumOfDirections_; ++cnt2)
        {
            // middle dots
            d_rib01 = nullptr;
            d_rib02 = nullptr;
            d_rib12 = nullptr;
            // search for ready ones
            for(size_t cntm=0; cntm!=DirTrigon2[cnt2][0]->midpointsNumber(); ++cntm )
            {
                if ( DirTrigon2[cnt2][0]->neighbor(cntm) == DirTrigon2[cnt2][1] )
                {
                    d_rib01=DirTrigon2[cnt2][0]->midpoint(cntm);
                    break;
                }
            }
            for(size_t cntm=0; cntm!=DirTrigon2[cnt2][0]->midpointsNumber(); ++cntm )
            {
                if ( DirTrigon2[cnt2][0]->neighbor(cntm) == DirTrigon2[cnt2][2] )
                {
                    d_rib02=DirTrigon2[cnt2][0]->midpoint(cntm);
                    break;
                }
            }
            for(size_t cntm=0; cntm!=DirTrigon2[cnt2][1]->midpointsNumber(); ++cntm )
            {
                if ( DirTrigon2[cnt2][1]->neighbor(cntm) == DirTrigon2[cnt2][2] )
                {
                    d_rib12=DirTrigon2[cnt2][1]->midpoint(cntm);
                    break;
                }
            }
            // new ones adding
            if (d_rib01 == nullptr)
            {
                DirDots[ NumOfReadyDirDots ].setMiddle(DirTrigon2[cnt2][0], DirTrigon2[cnt2][1] );
                d_rib01 = &(DirDots[ NumOfReadyDirDots ]);
                ++NumOfReadyDirDots;
                DirTrigon2[cnt2][0]->addMidpoint( DirTrigon2[cnt2][1], d_rib01 );
                DirTrigon2[cnt2][1]->addMidpoint( DirTrigon2[cnt2][0], d_rib01 );
            }
            if (d_rib02 == nullptr)
            {
                DirDots[ NumOfReadyDirDots ].setMiddle(DirTrigon2[cnt2][0], DirTrigon2[cnt2][2] );
                d_rib02 = &(DirDots[ NumOfReadyDirDots ]);
                ++NumOfReadyDirDots;
                DirTrigon2[cnt2][0]->addMidpoint( DirTrigon2[cnt2][2], d_rib02 );
                DirTrigon2[cnt2][2]->addMidpoint( DirTrigon2[cnt2][0], d_rib02 );
            }
            if (d_rib12 == nullptr)
            {
                DirDots[ NumOfReadyDirDots ].setMiddle(DirTrigon2[cnt2][1], DirTrigon2[cnt2][2] );
                d_rib12 = &(DirDots[ NumOfReadyDirDots ]);
                ++NumOfReadyDirDots;
                DirTrigon2[cnt2][1]->addMidpoint( DirTrigon2[cnt2][2], d_rib12 );
                DirTrigon2[cnt2][2]->addMidpoint( DirTrigon2[cnt2][1], d_rib12 );
            }
            // triangle splitting
            DirTrigon[ 4*cnt2   ].set(DirTrigon2[cnt2][0], d_rib01, d_rib02);
            DirTrigon[ 4*cnt2+1 ].set(DirTrigon2[cnt2][1], d_rib12, d_rib01);
            DirTrigon[ 4*cnt2+2 ].set(DirTrigon2[cnt2][2], d_rib02, d_rib12);
            DirTrigon[ 4*cnt2+3 ].set(d_rib01, d_rib12, d_rib02);
        }
        delete[] DirTrigon2;

        for (uint64_t cnt2=0; cnt2!=NumOfReadyDirDots; ++cnt2)
        {
            DirDots[cnt2].removeMidpoints();
        }

        NumOfDirections_ *= 4;
    }

    dots_ = new Vector3d[ NumOfDirections_ ];
    w_	  = new double[ NumOfDirections_ ];

    for (uint64_t cnt=0; cnt<NumOfDirections_; ++cnt)
    {
        dots_[cnt] = DirTrigon[cnt].median();
        w_[cnt] = DirTrigon[cnt].square()*NumOfDirections_ / 4. / PI;
    }
    delete[] DirDots;
    delete[] DirTrigon;
}
