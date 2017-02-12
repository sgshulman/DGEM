#ifndef _DIRECTIONS_HPP_
#define _DIRECTIONS_HPP_

// точка в пространстве
class DOT_N ;

// точка между двумя другими
class MID_DOT
{
	public:
		MID_DOT( void ) : d2_(nullptr), md_(nullptr) {}
		void Set( DOT_N* d2, DOT_N *md )
		{	d2_ = d2;	md_ = md;	} 
		DOT_N* d2( void )
		{	return d2_;	}
		DOT_N* md( void )
		{	return md_;	}
	private:
		DOT_N *d2_;
		DOT_N *md_;
};

class DOT_N 
{
	public:
		DOT_N( double x=0.0, double y=0.0, double z=0.0 ) : x_(x), y_(y), z_(z), NumOfMids_(0) {}; 
		void Get( double &x, double &y, double &z ) const
		{ x=x_; y=y_; z=z_;}
		void Set( double x, double y, double z ) 
		{ x_=x; y_=y; z_=z;	}
		double x( void ) const
		{	return x_;	}
		double y( void ) const
		{	return y_;	}
		double z( void ) const
		{	return z_;	}
		void Set(DOT_N *d1, DOT_N *d2)
		{
			double r, r1, r2;
			double m;
			
			x_ = ( d1->x() + d2->x() )/2.0;
			y_ = ( d1->y() + d2->y() )/2.0;
			z_ = ( d1->z() + d2->z() )/2.0;
				
			r1 = sqrt( d1->x()*d1->x() + d1->y()*d1->y() + d1->z()*d1->z() );
			r2 = sqrt( d2->x()*d2->x() + d2->y()*d2->y() + d2->z()*d2->z() );
			// сохранение R
			r = sqrt( x_*x_ + y_*y_ + z_*z_ );		
			m = ( r1+r2 )/( 2.0*r );
			x_*=m; y_*=m; z_*=m;
			NumOfMids_ = 0;
		}
		void Set1(DOT_N *d1, DOT_N *d2)
		{
			double r, r1, r2;
			double m;
			
			x_ = ( 1.9*d1->x() + 1.1*d2->x() )/3.0;
			y_ = ( 1.9*d1->y() + 1.1*d2->y() )/3.0;
			z_ = ( 1.9*d1->z() + 1.1*d2->z() )/3.0;
				
			r1 = sqrt( d1->x()*d1->x() + d1->y()*d1->y() + d1->z()*d1->z() );
			r2 = sqrt( d2->x()*d2->x() + d2->y()*d2->y() + d2->z()*d2->z() );
			// сохранение R
			r = sqrt( x_*x_ + y_*y_ + z_*z_ );		
			m = ( 2*r1+r2 )/( 3.0*r );
			x_*=m; y_*=m; z_*=m;
			NumOfMids_ = 0;
		}
		void Set2(DOT_N *d1, DOT_N *d2)
		{
			double r, r1, r2;
			double m;
			
			x_ = ( 1.1*d1->x() + 1.9*d2->x() )/3.0;
			y_ = ( 1.1*d1->y() + 1.9*d2->y() )/3.0;
			z_ = ( 1.1*d1->z() + 1.9*d2->z() )/3.0;
				
			r1 = sqrt( d1->x()*d1->x() + d1->y()*d1->y() + d1->z()*d1->z() );
			r2 = sqrt( d2->x()*d2->x() + d2->y()*d2->y() + d2->z()*d2->z() );
			// сохранение R
			r = sqrt( x_*x_ + y_*y_ + z_*z_ );		
			m = ( r1+2*r2 )/( 3.0*r );
			x_*=m; y_*=m; z_*=m;
			NumOfMids_ = 0;
		}
		uint32_t NumOfMids( void )
		{	return NumOfMids_;	}
		MID_DOT md( int i )
		{	return md_[ i ];	}
		void SetMid(DOT_N* d2, DOT_N *md)
		{
			md_[NumOfMids_].Set( d2, md );
			++NumOfMids_;
		}
		void ZeroMid( void )
		{	NumOfMids_ = 0;	}
	private:
		double x_, y_, z_;
		uint32_t NumOfMids_;
		MID_DOT md_[6];	
};

class DOT
{
	public:
		DOT( double x=0.0, double y=0.0, double z=0.0 ) : x_(x), y_(y), z_(z) {}; 
		DOT( DOT_N d ) : x_(d.x()), y_(d.y()), z_(d.z()) {};
		void Get( double &x, double &y, double &z ) 
		{ x=x_; y=y_; z=z_;}
		double x( void ) const
		{	return x_;	}
		double y( void ) const 
		{	return y_;	}
		double z( void ) const
		{	return z_;	}
	private:
		double x_, y_, z_;
};	

class DIRECTIONS {
	public:
		DIRECTIONS( uint32_t NumOfDirectionsLevels );
		~DIRECTIONS()
		{
			if ( Dots_ != nullptr )
			{
				delete[] Dots_;
				delete[] w_;
			}	
		}
		uint64_t NumOfDirections() const
		{	return ( NumOfDirections_ );	}
		void GetDirection( int index, double& x, double &y, double &z ) const
		{
			x = Dots_[ index ].x();
			y = Dots_[ index ].y();
			z = Dots_[ index ].z();
		}	
		double W( int index ) const
		{	return w_[ index ];	}
	private:
		DOT	*Dots_;	
		double *w_;
		uint64_t NumOfDirections_;
		
		DIRECTIONS ( DIRECTIONS const &);
		DIRECTIONS & operator =( DIRECTIONS const &);
};

#endif
