#include <math.h>
#include "grid.hpp"
#include "observers.hpp"
#include "model.hpp"
#include "photons.hpp"
#include "directions.hpp"

void Dustmat( double &p1, double &p2, double &p3, double &p4,
            double cost, double cost2, double pl, double pc,
            double sc, double hgg, double g2 )
{
// this calculates the elements of the phase matrix for a
//	simple representation of the mrn dust mixture using the algorithms
//	for the ultraviolet region due to richard l. white ap.j. 229, 954, 	1979.
    double phi, f, f2, c;

    p1 = (1 - g2)/pow( (1+g2-2*hgg*cost), 1.5);
    p2 = -pl*p1*(1-cost2)/(1+cost2);
    p3 = p1*2*cost/(1+cost2);

    if(cost >= 1.0) cost= 1.0;
    if(cost <=-1.0) cost=-1.0;

//     angle in degrees!
    phi=acos(cost)*180./3.1415926;
    f=3.13*phi*exp(-7.0*phi/180.);
//     now convert to radians
    f2=(phi+sc*f)*3.1415926/180.0;

    c=(cos(f2))*(cos(f2));
    p4 = -pc*p1*(1-c)/(1+c);
}

Photon::Photon( Position const &pos, double weight, int nscat )
{
    // Set position
    pos_=pos;
    // Set direction
    double v, u;

     v = ran.Get();
     u = ran.Get();

    dir_ = Direction( 2.0*PI*u, acos( 2*v-1.0 ) );
    // Set number of sctterings and weight
    nscat_ = nscat;
    weight_= weight;
    // Set Stokes fluxes
    fi_=1.0;
    fq_=0.0;
    fu_=0.0;
    fv_=0.0;
}

Photon::Photon( Position const &pos, Direction const &dir, double weight, int nscat, double fi, double fq, double fu, double fv )
{
    pos_ = pos;
    dir_ = dir;
    nscat_ = nscat;
    weight_ = weight;
    // Set Stokes fluxes
    fi_=fi;
    fq_=fq;
    fu_=fu;
    fv_=fv;
}

double Photon::Scatt( Model const &m, Direction const & dir )
{
    double calpha; // cos(alpha), where alpha is angle between incident
                    // and outgoing (i.e., observed) photon direction
    double hgfrac;
    calpha=dir_.nx()*dir.nx()+dir_.ny()*dir.ny()+dir_.nz()*dir.nz();

    //weight photon for isotropic, Thomson, or HG scattering
    // hgfrac=0.5*(1.0+calpha*calpha)  ;       // Thomson
    // hgfrac=1./4/3.1415926    ;              // isotropic
     hgfrac=(1.0-m.g2())/pow((1.0+m.g2()-2.*m.hgg()*calpha),1.5)/4/3.1415926; // HG
       
     Stokes( m, dir, calpha, true );
 //    dir_ = dir;
     return hgfrac;
}

void Photon::Scatt( Model const &m, Directions const &dirs, Grid const &grid, std::vector<Observer>& observers)
{
    if (nscat_ == m.MonteCarloStart() )
    {
        Photon ph(*this);
        int tflag = 0;
        ph.Stokes( m, Direction(), 0.0, false );
        tflag = grid.TauInt2( ph );
        ph.nscat()+=1;
        while ( !tflag && ( ph.nscat() <= m.nscat() ) )
        {
            ph.weight() *= m.albedo();
            // Do peeling off and project weighted photons into image
            // учитыается нерассеяшийся свет от каждой точки рассеяния и последующие рассеяния, пока фотон не изыдет
            for (Observer& observer : observers)
            {
                grid.Peeloff(ph, observer, m);
            }

            // Scatter photon into new direction and update Stokes parameters
            ph.Stokes( m, Direction(), 0.0, false );
            ph.nscat()+=1;
            if (ph.nscat() > m.nscat()) break;
            // Find next scattering location
            tflag = grid.TauInt2( ph );
        }
    } else {
        double calpha;
        double hgfrac;
        // normalizing
        double sum=0;
        for (size_t j=0; j!=dirs.NumOfDirections(); ++j)
        {
            double x, y, z;
            dirs.GetDirection( j, x, y, z );
            calpha = dir_.nx()*x+dir_.ny()*y+dir_.nz()*z;
            sum += dirs.W( j )*(1.0-m.g2())/pow((1.0+m.g2()-2.*m.hgg()*calpha),1.5);
        }

        // randomized grid experiment
        /*double v = ran.Get();
        double u = ran.Get();
        double alpha =  2.0*PI*u;
        double beta = acos( 2*v-1.0 );*/
        // randomized grid experiment

        for (size_t j=0; j!=dirs.NumOfDirections(); ++j)
        {
            // Release photon from point source
            double x, y, z;
            dirs.GetDirection( j, x, y, z );

            // randomized grid experiment
            /*double t1, t2;
            t1 = x*cos(alpha)+y*sin(alpha);
            t2 =-x*sin(alpha)+y*cos(alpha);
            x = t1;
            y = t2;
            t1 = x*cos(beta)+z*sin(beta);
            t2 =-x*sin(beta)+z*cos(beta);
            x = t1;
            z = t2;*/
            // randomized grid experiment

            calpha = dir_.nx()*x+dir_.ny()*y+dir_.nz()*z;
            hgfrac=(1.0-m.g2())/pow((1.0+m.g2()-2.*m.hgg()*calpha),1.5);
            Photon ph0(pos_, dir_, weight_*dirs.W( j )*hgfrac/sum, nscat_+1, fi_, fq_, fu_, fv_ );
            ph0.Stokes(m, Direction(x, y, z), calpha, true);

            // Find optical depth, tau1, to edge of grid
            //double tau1 = grid.TauFind( ph0 );
            //if ( tau1 < m.taumin() ) continue;

            double w = 1.0 / m.NumOfSecondaryScatterings() ;
            Position spos = pos_;
            double tauold = 0.0, tau = 0.0;

            // Loop over scattering dots
            for (size_t s=0; s!=m.NumOfSecondaryScatterings(); ++s)
            {
                Photon ph( spos, ph0.dir(), ph0.weight()*w, nscat_+1, ph0.fi(), ph0.fq(), ph0.fu(), ph0.fv() );
                // Force photon to scatter at optical depth tau before edge of grid
                tauold = tau;
                tau=-log( 1.0-0.5*w*(2*s+1) );
                // Find scattering location of tau
                if( grid.TauInt( ph, tau, tauold ) )
                {
                    break;
                } else {
                    spos = ph.pos();

                    // Photon scattering
                    ph.weight() *= m.albedo();

                    for (Observer& observer : observers)
                    {
                        grid.Peeloff(ph, observer, m);
                    }

                    if (ph.nscat() < m.nscat() ) ph.Scatt( m, dirs, grid, observers );
                }
            }
        }
    }
}

// Stokes vector changes
// spherical trigonometry is used
void Photon::Stokes( Model const &m, Direction const &dir, double calpha, bool fDir )
{
    double a11,a12,a13,a21,a22,a23,a24,a31,a32,a33,a34;
        double a42,a43,a44;
    double phi=0.0, cost, sint;

    double wght = fi_;
    fi_ /= wght;
    fq_ /= wght;
    fu_ /= wght;
    fv_ /= wght;

    double costp = dir_.cost();
    double sintp = dir_.sint();
    double phip = dir_.phi();

    // dust scattering
    double cosTh, cos2Th; //bmu, bmu2;
    if (fDir)
    {
        cosTh = calpha;
    } else {
        cosTh=((1.0+m.g2())-pow( (( 1.0-m.g2() )/( 1.0-m.hgg()+2.0*m.hgg()*ran.Get())), 2 ))/( 2.0*m.hgg() );
    }
    if ( cosTh > 1.0 )
    {
        cosTh=1.0;
    }else if( cosTh < -1.0 ) {
        cosTh=-1.0;
    }
    cos2Th = cosTh*cosTh;

    double p1, p2, p3, p4;
    Dustmat(p1,p2,p3,p4,cosTh,cos2Th,m.pl(),m.pc(),m.sc(),m.hgg(),m.g2());
    double a = p1;

    double sinTh = sqrt(1.0-cos2Th);//sinbt = sqrt(1.0-bmu2);

    double ri1;
    if (fDir)
    {
        double ophi = dir.phi();
        if (ophi > 3.1415926) ophi -= 2*3.1415926;
        double cosi1 = (dir.cost()-dir_.cost()*cosTh)/(dir_.sint()*sinTh);
        double sini1 = sin(dir_.phi()-ophi-3.1415926)*dir.sint()/sinTh;
        ri1 = atan2(sini1, cosi1)+3.1415926;
    } else {
        ri1=2*3.1415926*ran.Get();
    }


    if(ri1 >= 3.1415926)
    {
        double ri3 = 3.1415926*2-ri1;
        double cosi3=cos(ri3);
        double sini3=sin(ri3);
        double sin2i3=2.0*sini3*cosi3;
        double cos2i3=2.0*cosi3*cosi3-1.0;
        a11=p1;
        a12=p2*cos2i3;
        a13=p2*sin2i3;

        if(cosTh == 1.0) return;
        if(cosTh ==-1.0)
        {
            fu_=-fu_;
            return;
        }
        if (fDir)
        {
            cost = dir.cost();
            sint = dir.sint();
        } else {
            cost = costp*cosTh+sintp*sinTh*cosi3;
            sint = fabs(sqrt(1.0-cost*cost));
        }
        double sini2=0.0, cosi2=1.0;
        if(cost < 1.0 && cost > -1.0)
        {
            sini2=sini3*sintp/sint;
            double bott=sint*sinTh;
            cosi2=costp/bott-cost*cosTh/bott;
        } else if(cost >= 1.0) {
            cosi2=-1.0;
        }

        if (!fDir)
        {
             double cosdph=-cosi2*cosi3+sini2*sini3*cosTh;
             if(cosdph > 1.0) cosdph=1.0;
             if(cosdph <-1.0) cosdph=-1.0;

             phi=phip+acos(cosdph);
             if(phi > 2*3.1415926 ) phi=phi-2*3.1415926;
             if(phi < 0.0)          phi=phi+2*3.1415926;
        }
        double sin2i2=2.0*sini2*cosi2;
        double cos2i2=2.0*cosi2*cosi2-1.0;
        double sin2=sin2i2*sin2i3;
        double cos2=cos2i2*cos2i3;
        double sin2cos1=sin2i2*cos2i3;
        double cos2sin1=cos2i2*sin2i3;

         a21=p2*cos2i2;
         a22=p1*cos2-p3*sin2;
         a23=p1*cos2sin1+p3*sin2cos1;
         a24=-p4*sin2i2;
         a31=-p2*sin2i2;
         a32=-p1*sin2cos1-p3*cos2sin1;
         a33=-p1*sin2+p3*cos2;
         a34=-p4*cos2i2;
         a42=-p4*sin2i3;
         a43=p4*cos2i3;
         a44=p3;
    } else {
        double cosi1=cos(ri1);
        double sini1=sin(ri1);
        double sin2i1=2.0*sini1*cosi1;
        double cos2i1=2.0*cosi1*cosi1-1.0;
        a11=p1;
        a12=p2*cos2i1;
        a13=-p2*sin2i1;

        if(cosTh == 1.0) return;
        if(cosTh ==-1.0)
        {
            fu_=-fu_;
            return;
        }
        if ( fDir )
        {
            cost = dir.cost();
            sint = dir.sint();
        } else {
            cost=costp*cosTh+sintp*sinTh*cosi1;
            sint=fabs(sqrt(1.0-cost*cost));
        }
        double sini2=0.0, cosi2=1.0;
        if(cost < 1.0 && cost > -1.0)
         {
            sini2=sini1*sintp/sint;
            double bott=sint*sinTh;
            cosi2=costp/bott-cost*cosTh/bott;
         } else if(cost >= 1.0) {
            cosi2=-1.0;
        }

        if( !fDir )
        {
            double cosdph=-cosi1*cosi2+sini1*sini2*cosTh;
            if( cosdph > 1.0 ) cosdph= 1.0;
            if( cosdph <-1.0 ) cosdph=-1.0;

            phi=phip-acos(cosdph);
             if(phi > 2*3.1415926 ) phi=phi-2*3.1415926;
             if(phi < 0.0)          phi=phi+2*3.1415926;
        }
        double sin2i2=2.0*sini2*cosi2;
        double cos2i2=2.0*cosi2*cosi2-1.0;
        double sin2=sin2i2*sin2i1;
        double cos2=cos2i2*cos2i1;
        double sin2cos1=sin2i2*cos2i1;
        double cos2sin1=cos2i2*sin2i1;

         a21=p2*cos2i2;
         a22=p1*cos2-p3*sin2;
         a23=-p1*cos2sin1-p3*sin2cos1;
         a24=p4*sin2i2;
         a31=p2*sin2i2;
         a32=p1*sin2cos1+p3*cos2sin1;
         a33=-p1*sin2+p3*cos2;
         a34=-p4*cos2i2;
         a42=p4*sin2i1;
         a43=p4*cos2i1;
         a44=p3;
    }

    double si=(a11*fi_+a12*fq_+a13*fu_)/a;
    double sq=(a21*fi_+a22*fq_+a23*fu_+a24*fv_)/a;
    double su=(a31*fi_+a32*fq_+a33*fu_+a34*fv_)/a;
    double sv=(a42*fq_+a43*fu_+a44*fv_)/a;

    fi_=si*wght;
    fq_=sq*wght;
    fu_=su*wght;
    fv_=sv*wght;

    if( fDir )
    {
        dir_ = dir;
    } else {
        double cosp=cos(phi);
        double sinp=sin(phi);
        dir_ = Direction(sint*cosp, sint*sinp, cost);
    }
}


