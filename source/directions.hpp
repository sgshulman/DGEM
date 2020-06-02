#ifndef DIRECTIONS_HPP_
#define DIRECTIONS_HPP_

// точка в пространстве
class DotN ;

// точка между двумя другими
class MiddleDot
{
	public:
		MiddleDot( void ) : d2_(nullptr), md_(nullptr) {}
		void Set( DotN* d2, DotN *md )
		{	d2_ = d2;	md_ = md;	} 
		DotN* d2( void )
		{	return d2_;	}
		DotN* md( void )
		{	return md_;	}
	private:
		DotN *d2_;
		DotN *md_;
};

class DotN 
{
	public:
		DotN( double x=0.0, double y=0.0, double z=0.0 ) : x_(x), y_(y), z_(z), NumOfMids_(0) {}; 
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
		void Set(DotN *d1, DotN *d2)
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
		void Set1(DotN *d1, DotN *d2)
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
		void Set2(DotN *d1, DotN *d2)
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
		MiddleDot md( int i )
		{	return md_[ i ];	}
		void SetMid(DotN* d2, DotN *md)
		{
			md_[NumOfMids_].Set( d2, md );
			++NumOfMids_;
		}
		void ZeroMid( void )
		{	NumOfMids_ = 0;	}
	private:
		double x_, y_, z_;
		uint32_t NumOfMids_;
		MiddleDot md_[6];	
};

class Dot
{
	public:
		Dot( double x=0.0, double y=0.0, double z=0.0 ) : x_(x), y_(y), z_(z) {}; 
		Dot( DotN d ) : x_(d.x()), y_(d.y()), z_(d.z()) {};
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

class Directions {
	public:
		Directions( uint32_t NumOfDirectionsLevels );
		~Directions()
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
		Dot	*Dots_;	
		double *w_;
		uint64_t NumOfDirections_;
		
		Directions ( Directions const &);
		Directions & operator =( Directions const &);
};

#endif
