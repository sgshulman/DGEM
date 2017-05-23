#include <iostream>
#include <cmath>
#include <cstdio>

int main( int argc, char *argv[] )
{
	if( argc!=3 ) return 0;
	int Nx_ = 200;
	int Ny_ = 200;
	
	FILE *f1=NULL;
	FILE *f2=NULL;
	FILE *fd1=NULL;
	FILE *fd2=NULL;
	
	f1= fopen(argv[1], "r");
	f2= fopen(argv[2], "r");
	fd1= fopen("dif_abs.dat", "w");
	fd2= fopen("dif_rel.dat", "w");
	double sum=0.0, sum2=0.0, totsum=0;
	double a, b, tmp;
	int num=0;
	for (long int x=0; x!=Nx_; ++x)
	{
		for (long int y=0; y!=Ny_; ++y)
		{
			fscanf(f1, "%lg", &a);
			fscanf(f2, "%lg", &b);
			tmp = fabs(a-b);	
			fprintf(fd1, "%lg\t", tmp);
			if (0.5*fabs(a+b) < 0.07)
			{
				sum+=tmp;
				totsum += 0.5*fabs(a+b);
			} else {
			    printf("ZERO\n");
			}
			if (a!=0 && b!=0)
			{
				++num;
				sum2+=2*tmp/fabs(a+b);
				fprintf(fd2, "%lg\t", 2*tmp/fabs(a+b));
			} else {
				fprintf(fd2, "0\t");
			}
			
		}
		fprintf(fd1, "\n");
		fprintf(fd2, "\n");
	}
	printf("sum = %lg sum/n=%lg sum2/n=%lg n =%i \n refsum = %lg\n", sum, sum/num, sum2/num, num, sum/totsum );
	fclose( f1 );
	fclose( f2 );
	fclose( fd1 );
	fclose( fd2 );
}
