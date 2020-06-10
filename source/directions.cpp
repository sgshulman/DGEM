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
            PointWithNeighbors(double x = 0.0, double y = 0.0, double z = 0.0) : x_(x), y_(y), z_(z), midpointsNumber_(0)
            {};

            void Get(double &x, double &y, double &z) const
            {
                x = x_;
                y = y_;
                z = z_;
            }

            void Set(double x, double y, double z)
            {
                x_ = x;
                y_ = y;
                z_ = z;
            }

            double x(void) const
            { return x_; }

            double y(void) const
            { return y_; }

            double z(void) const
            { return z_; }

            void Set(PointWithNeighbors const *d1, PointWithNeighbors const *d2)
            {
                double r, r1, r2;
                double m;

                x_ = (d1->x() + d2->x()) / 2.0;
                y_ = (d1->y() + d2->y()) / 2.0;
                z_ = (d1->z() + d2->z()) / 2.0;

                r1 = sqrt(d1->x() * d1->x() + d1->y() * d1->y() + d1->z() * d1->z());
                r2 = sqrt(d2->x() * d2->x() + d2->y() * d2->y() + d2->z() * d2->z());
                // сохранение R
                r = sqrt(x_ * x_ + y_ * y_ + z_ * z_);
                m = (r1 + r2) / (2.0 * r);
                x_ *= m;
                y_ *= m;
                z_ *= m;
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
            double x_, y_, z_;
            uint32_t midpointsNumber_;
            PointWithNeighbors const *neighbors_[6];
            PointWithNeighbors *midpoints_[6];
    };


    // Element of a mesh on the union sphere
    class SphericalTriangle
    {
        public:
            SphericalTriangle()
            {}

            void set(PointWithNeighbors *dot1, PointWithNeighbors *dot2, PointWithNeighbors *dot3)
            {
                dots_[0] = dot1;
                dots_[1] = dot2;
                dots_[2] = dot3;
            }

            Vector3d median()
            {
                double x, y, z, r;
                x = (dots_[0]->x() + dots_[1]->x() + dots_[2]->x()) / 3.0;
                y = (dots_[0]->y() + dots_[1]->y() + dots_[2]->y()) / 3.0;
                z = (dots_[0]->z() + dots_[1]->z() + dots_[2]->z()) / 3.0;
                r = sqrt(x * x + y * y + z * z);
                x /= r;
                y /= r;
                z /= r;
                return {x, y, z};
            }

            double square(void)
            {
                double r1, r2, r3;

                r1 = sqrt(dots_[0]->x() * dots_[0]->x() + dots_[0]->y() * dots_[0]->y() + dots_[0]->z() * dots_[0]->z());
                r2 = sqrt(dots_[1]->x() * dots_[1]->x() + dots_[1]->y() * dots_[1]->y() + dots_[1]->z() * dots_[1]->z());
                r3 = sqrt(dots_[2]->x() * dots_[2]->x() + dots_[2]->y() * dots_[2]->y() + dots_[2]->z() * dots_[2]->z());

                double mix = dots_[0]->x() * dots_[1]->y() * dots_[2]->z() + dots_[0]->y() * dots_[1]->z() * dots_[2]->x() +
                             dots_[0]->z() * dots_[1]->x() * dots_[2]->y() - dots_[0]->z() * dots_[1]->y() * dots_[2]->x() -
                             dots_[0]->y() * dots_[1]->x() * dots_[2]->z() - dots_[0]->x() * dots_[1]->z() * dots_[2]->y();

                double r1r2 = dots_[0]->x() * dots_[1]->x() + dots_[0]->y() * dots_[1]->y() + dots_[0]->z() * dots_[1]->z();
                double r2r3 = dots_[1]->x() * dots_[2]->x() + dots_[1]->y() * dots_[2]->y() + dots_[1]->z() * dots_[2]->z();
                double r1r3 = dots_[0]->x() * dots_[2]->x() + dots_[0]->y() * dots_[2]->y() + dots_[0]->z() * dots_[2]->z();

                return 2 * atan(mix / (r1 * r2 * r3 + r1r2 * r3 + r2r3 * r1 + r1r3 * r2));
            }

            PointWithNeighbors *operator[](int i)
            { return dots_[i]; }

        private:
            PointWithNeighbors *dots_[3]{};
    };
}

// directions
Directions::Directions( uint32_t NumOfDirectionsLevels )
{
    NumOfDirections_ = 20;
    unsigned int NumOfReadyDirDots;
    int NumOfDirDots = 12;
    for(size_t cnt=1; cnt<NumOfDirectionsLevels; cnt++)
    {
        NumOfDirDots = NumOfDirDots+(NumOfDirections_*3)/2;
        NumOfDirections_ *= 4;
    }
    NumOfDirections_ = 20;

    PointWithNeighbors *DirDots = new PointWithNeighbors[ NumOfDirDots ];

    SphericalTriangle *DirTrigon = new SphericalTriangle[ NumOfDirections_ ];
    SphericalTriangle *DirTrigon2;

    PointWithNeighbors *d_rib01, *d_rib02, *d_rib12;
    // initial sphere conditions
    (DirDots[0]).Set(  0.0, 0.0, sqrt(5.0)/2.0 );
    for (size_t cnt =0; cnt<5; cnt++)
    {
        (DirDots[cnt+1]).Set(  cos(cnt*72*PI/180),  sin(cnt*72*PI/180),  0.5 );
        (DirDots[cnt+6]).Set(  cos((36+cnt*72)*PI/180),  sin((36+cnt*72)*PI/180),  -0.5 );
    }
    (DirDots[11]).Set( 0.0, 0.0, -sqrt(5.0)/2.0   );
    NumOfReadyDirDots = 12;
    // triangles
    (DirTrigon[0]).set( &(DirDots[0]), &(DirDots[1]), &(DirDots[2]) );
    (DirTrigon[1]).set( &(DirDots[0]), &(DirDots[2]), &(DirDots[3]) );
    (DirTrigon[2]).set( &(DirDots[0]), &(DirDots[3]), &(DirDots[4]) );
    (DirTrigon[3]).set( &(DirDots[0]), &(DirDots[4]), &(DirDots[5]) );
    (DirTrigon[4]).set( &(DirDots[0]), &(DirDots[5]), &(DirDots[1]) );
    (DirTrigon[5]).set( &(DirDots[10]), &(DirDots[11]), &(DirDots[6]) );
    (DirTrigon[6]).set( &(DirDots[11]), &(DirDots[7]), &(DirDots[6]) );
    (DirTrigon[7]).set( &(DirDots[11]), &(DirDots[8]), &(DirDots[7]) );
    (DirTrigon[8]).set( &(DirDots[11]), &(DirDots[9]), &(DirDots[8]) );
    (DirTrigon[9]).set( &(DirDots[11]), &(DirDots[10]), &(DirDots[9]) );
    (DirTrigon[10]).set( &(DirDots[6]), &(DirDots[1]), &(DirDots[10]) );
    (DirTrigon[11]).set( &(DirDots[2]), &(DirDots[6]), &(DirDots[7]) );
    (DirTrigon[12]).set( &(DirDots[3]), &(DirDots[7]), &(DirDots[8]) );
    (DirTrigon[13]).set( &(DirDots[4]), &(DirDots[8]), &(DirDots[9]) );
    (DirTrigon[14]).set( &(DirDots[5]), &(DirDots[9]), &(DirDots[10]) );
    (DirTrigon[15]).set( &(DirDots[6]), &(DirDots[2]), &(DirDots[1]) );
    (DirTrigon[16]).set( &(DirDots[7]), &(DirDots[3]), &(DirDots[2]) );
    (DirTrigon[17]).set( &(DirDots[8]), &(DirDots[4]), &(DirDots[3]) );
    (DirTrigon[18]).set( &(DirDots[9]), &(DirDots[5]), &(DirDots[4]) );
    (DirTrigon[19]).set( &(DirDots[10]), &(DirDots[1]), &(DirDots[5]) );

    for(size_t cnt=1; cnt<NumOfDirectionsLevels; cnt++)
    {
        // new triangles
        DirTrigon2 = DirTrigon;
        DirTrigon = new SphericalTriangle[ NumOfDirections_*4 ];
        // splitting
        for(size_t cnt2 = 0; cnt2 < NumOfDirections_; cnt2++)
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
                DirDots[ NumOfReadyDirDots ].Set(DirTrigon2[cnt2][0], DirTrigon2[cnt2][1] );
                d_rib01 = &(DirDots[ NumOfReadyDirDots ]);
                NumOfReadyDirDots++;
                DirTrigon2[cnt2][0]->addMidpoint( DirTrigon2[cnt2][1], d_rib01 );
                DirTrigon2[cnt2][1]->addMidpoint( DirTrigon2[cnt2][0], d_rib01 );
            }
            if (d_rib02 == nullptr)
            {
                DirDots[ NumOfReadyDirDots ].Set(DirTrigon2[cnt2][0], DirTrigon2[cnt2][2] );
                d_rib02 = &(DirDots[ NumOfReadyDirDots ]);
                NumOfReadyDirDots++;
                DirTrigon2[cnt2][0]->addMidpoint( DirTrigon2[cnt2][2], d_rib02 );
                DirTrigon2[cnt2][2]->addMidpoint( DirTrigon2[cnt2][0], d_rib02 );
            }
            if (d_rib12 == nullptr)
            {
                DirDots[ NumOfReadyDirDots ].Set(DirTrigon2[cnt2][1], DirTrigon2[cnt2][2] );
                d_rib12 = &(DirDots[ NumOfReadyDirDots ]);
                NumOfReadyDirDots++;
                DirTrigon2[cnt2][1]->addMidpoint( DirTrigon2[cnt2][2], d_rib12 );
                DirTrigon2[cnt2][2]->addMidpoint( DirTrigon2[cnt2][1], d_rib12 );
            }
            // triangle splitting
            DirTrigon[ 4*cnt2   ].set(DirTrigon2[cnt2][0], d_rib01, d_rib02);
            DirTrigon[ 4*cnt2+1 ].set(DirTrigon2[cnt2][1], d_rib12, d_rib01);
            DirTrigon[ 4*cnt2+2 ].set(DirTrigon2[cnt2][2], d_rib02, d_rib12);
            DirTrigon[ 4*cnt2+3 ].set( d_rib01, d_rib12, d_rib02 );
        }
        delete[] DirTrigon2;
        DirTrigon2 = nullptr;

        for (unsigned int cnt2=0; cnt2!=NumOfReadyDirDots; ++cnt2)
            DirDots[ cnt2 ].removeMidpoints();

        NumOfDirections_ *= 4;
    }

    dots_ = new Vector3d[ NumOfDirections_ ];
    w_	  = new double[ NumOfDirections_ ];

    for (size_t cnt=0; cnt<NumOfDirections_; cnt++)
    {
        dots_[cnt] = DirTrigon[cnt].median();
        w_[cnt] = DirTrigon[cnt].square()*NumOfDirections_/ 4. / PI;
    }
    delete[] DirDots;
    delete[] DirTrigon;
}
