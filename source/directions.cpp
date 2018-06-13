#include <cmath>
#include <cstdio>
#include "model.hpp"
#include "directions.hpp"

// triangle
class MATTER 
{
	public:
		MATTER()
		{
			for (size_t cnt = 0; cnt<3; cnt++)
				dots_[ cnt ] = nullptr;
		}
		// init
		void Set(DOT_N *dot1, DOT_N *dot2, DOT_N *dot3 )
		{	
			dots_[0] = dot1;
			dots_[1] = dot2; 
			dots_[2] = dot3;
		}
		DOT_N Median( void )
		{
			double x, y, z, r;
			x = ( dots_[0]->x() + dots_[1]->x() + dots_[2]->x() )/3.0;
			y = ( dots_[0]->y() + dots_[1]->y() + dots_[2]->y() )/3.0;
			z = ( dots_[0]->z() + dots_[1]->z() + dots_[2]->z() )/3.0;
			r = sqrt( x*x + y*y + z*z );
			x /= r;
			y /= r;
			z /= r;
			return DOT_N(x, y, z);
		}
		double Square( void )
		{
			double r1, r2, r3;
			
			r1 = sqrt( dots_[0]->x()*dots_[0]->x()+dots_[0]->y()*dots_[0]->y()+dots_[0]->z()*dots_[0]->z() );
			r2 = sqrt( dots_[1]->x()*dots_[1]->x()+dots_[1]->y()*dots_[1]->y()+dots_[1]->z()*dots_[1]->z() );
			r3 = sqrt( dots_[2]->x()*dots_[2]->x()+dots_[2]->y()*dots_[2]->y()+dots_[2]->z()*dots_[2]->z() );
			
			double mix = dots_[0]->x()*dots_[1]->y()*dots_[2]->z() + dots_[0]->y()*dots_[1]->z()*dots_[2]->x() + 
					dots_[0]->z()*dots_[1]->x()*dots_[2]->y() - dots_[0]->z()*dots_[1]->y()*dots_[2]->x() - 
					dots_[0]->y()*dots_[1]->x()*dots_[2]->z() - dots_[0]->x()*dots_[1]->z()*dots_[2]->y();
					
			double r1r2 = dots_[0]->x()*dots_[1]->x()+dots_[0]->y()*dots_[1]->y()+dots_[0]->z()*dots_[1]->z();
			double r2r3 = dots_[1]->x()*dots_[2]->x()+dots_[1]->y()*dots_[2]->y()+dots_[1]->z()*dots_[2]->z();
			double r1r3 = dots_[0]->x()*dots_[2]->x()+dots_[0]->y()*dots_[2]->y()+dots_[0]->z()*dots_[2]->z();
			
			return 2*atan( mix/(r1*r2*r3+r1r2*r3+r2r3*r1+r1r3*r2) );
		}
		DOT_N* & operator []( int i )
		{	return dots_[ i ];	}
	private:
		DOT_N*	dots_[3];
};

// directions
DIRECTIONS::DIRECTIONS( uint32_t NumOfDirectionsLevels )
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
				
	DOT_N *DirDots = new DOT_N[ NumOfDirDots ];
	
	MATTER *DirTrigon = new MATTER[ NumOfDirections_ ];
	MATTER *DirTrigon2;
			
	DOT_N *d_rib01, *d_rib02, *d_rib12;
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
	(DirTrigon[0]).Set( &(DirDots[0]), &(DirDots[1]), &(DirDots[2]) );
	(DirTrigon[1]).Set( &(DirDots[0]), &(DirDots[2]), &(DirDots[3]) );
	(DirTrigon[2]).Set( &(DirDots[0]), &(DirDots[3]), &(DirDots[4]) );
	(DirTrigon[3]).Set( &(DirDots[0]), &(DirDots[4]), &(DirDots[5]) );
	(DirTrigon[4]).Set( &(DirDots[0]), &(DirDots[5]), &(DirDots[1]) );
	(DirTrigon[5]).Set( &(DirDots[10]), &(DirDots[11]), &(DirDots[6]) );
	(DirTrigon[6]).Set( &(DirDots[11]), &(DirDots[7]), &(DirDots[6]) );
	(DirTrigon[7]).Set( &(DirDots[11]), &(DirDots[8]), &(DirDots[7]) );
	(DirTrigon[8]).Set( &(DirDots[11]), &(DirDots[9]), &(DirDots[8]) );
	(DirTrigon[9]).Set( &(DirDots[11]), &(DirDots[10]), &(DirDots[9]) );
	(DirTrigon[10]).Set( &(DirDots[6]), &(DirDots[1]), &(DirDots[10]) );
	(DirTrigon[11]).Set( &(DirDots[2]), &(DirDots[6]), &(DirDots[7]) );
	(DirTrigon[12]).Set( &(DirDots[3]), &(DirDots[7]), &(DirDots[8]) );
	(DirTrigon[13]).Set( &(DirDots[4]), &(DirDots[8]), &(DirDots[9]) );
	(DirTrigon[14]).Set( &(DirDots[5]), &(DirDots[9]), &(DirDots[10]) );
	(DirTrigon[15]).Set( &(DirDots[6]), &(DirDots[2]), &(DirDots[1]) );
	(DirTrigon[16]).Set( &(DirDots[7]), &(DirDots[3]), &(DirDots[2]) );
	(DirTrigon[17]).Set( &(DirDots[8]), &(DirDots[4]), &(DirDots[3]) );
	(DirTrigon[18]).Set( &(DirDots[9]), &(DirDots[5]), &(DirDots[4]) );
	(DirTrigon[19]).Set( &(DirDots[10]), &(DirDots[1]), &(DirDots[5]) );
			
	for(size_t cnt=1; cnt<NumOfDirectionsLevels; cnt++)
	{		
		// new triangles
		DirTrigon2 = DirTrigon;
		DirTrigon = new MATTER[ NumOfDirections_*4 ];
		// splitting
		for(size_t cnt2 = 0; cnt2 < NumOfDirections_; cnt2++)
		{
			// midle dots
			d_rib01 = nullptr;
			d_rib02 = nullptr;
			d_rib12 = nullptr;
			// search for ready ones
			for(size_t cntm=0; cntm!=DirTrigon2[cnt2][0]->NumOfMids(); ++cntm )
			{
				if ( DirTrigon2[cnt2][0]->md(cntm).d2() == DirTrigon2[cnt2][1] )
				{
					d_rib01=DirTrigon2[cnt2][0]->md(cntm).md();
					break;
				}
			}
			for(size_t cntm=0; cntm!=DirTrigon2[cnt2][0]->NumOfMids(); ++cntm )
			{
				if ( DirTrigon2[cnt2][0]->md(cntm).d2() == DirTrigon2[cnt2][2] )
				{
					d_rib02=DirTrigon2[cnt2][0]->md(cntm).md();
					break;
				}
			}
			for(size_t cntm=0; cntm!=DirTrigon2[cnt2][1]->NumOfMids(); ++cntm )
			{
				if ( DirTrigon2[cnt2][1]->md(cntm).d2() == DirTrigon2[cnt2][2] )
				{
					d_rib12=DirTrigon2[cnt2][1]->md(cntm).md();
					break;
				}
			}
			// new ones adding
			if (d_rib01 == nullptr)
			{
				DirDots[ NumOfReadyDirDots ].Set(DirTrigon2[cnt2][0], DirTrigon2[cnt2][1] );
				d_rib01 = &(DirDots[ NumOfReadyDirDots ]);
				NumOfReadyDirDots++;
				DirTrigon2[cnt2][0]->SetMid( DirTrigon2[cnt2][1], d_rib01 );
				DirTrigon2[cnt2][1]->SetMid( DirTrigon2[cnt2][0], d_rib01 );
			}
			if (d_rib02 == nullptr)
			{
				DirDots[ NumOfReadyDirDots ].Set(DirTrigon2[cnt2][0], DirTrigon2[cnt2][2] );
				d_rib02 = &(DirDots[ NumOfReadyDirDots ]);
				NumOfReadyDirDots++;
				DirTrigon2[cnt2][0]->SetMid( DirTrigon2[cnt2][2], d_rib02 );
				DirTrigon2[cnt2][2]->SetMid( DirTrigon2[cnt2][0], d_rib02 );
			}
			if (d_rib12 == nullptr)
			{
				DirDots[ NumOfReadyDirDots ].Set(DirTrigon2[cnt2][1], DirTrigon2[cnt2][2] );
				d_rib12 = &(DirDots[ NumOfReadyDirDots ]);
				NumOfReadyDirDots++;
				DirTrigon2[cnt2][1]->SetMid( DirTrigon2[cnt2][2], d_rib12 );
				DirTrigon2[cnt2][2]->SetMid( DirTrigon2[cnt2][1], d_rib12 );
			}
			// triangle splitting
			DirTrigon[ 4*cnt2   ].Set(DirTrigon2[cnt2][0], d_rib01, d_rib02);
			DirTrigon[ 4*cnt2+1 ].Set(DirTrigon2[cnt2][1], d_rib12, d_rib01);
			DirTrigon[ 4*cnt2+2 ].Set(DirTrigon2[cnt2][2], d_rib02, d_rib12);
			DirTrigon[ 4*cnt2+3 ].Set( d_rib01, d_rib12, d_rib02 );			
		}
		delete[] DirTrigon2;
		DirTrigon2 = nullptr;

		for (unsigned int cnt2=0; cnt2!=NumOfReadyDirDots; ++cnt2)
			DirDots[ cnt2 ].ZeroMid();
			
		NumOfDirections_ *= 4;
	}				
		
	Dots_ = new DOT[ NumOfDirections_ ];
	w_	  = new double[ NumOfDirections_ ];

	for (size_t cnt=0; cnt<NumOfDirections_; cnt++)
	{	
		Dots_[cnt]=DOT(DirTrigon[cnt].Median());
		w_[cnt] = DirTrigon[cnt].Square()*NumOfDirections_/4/PI;
	}
	delete[] DirDots;		
	delete[] DirTrigon;
}
