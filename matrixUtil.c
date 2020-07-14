#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrixUtil.h"

#if 0
int main() 
{
	double *m;
	m = ( double * ) malloc( 100 * sizeof( double ));
	
	readMatrixFromFile( "test_cases/x10", m );

	int i = 0;
//	for( i = 0; i < 100; i ++ )
//		printf("%lf	", *( m + i ));

	free( m );
	return 0;
}
#endif

void readMatrixFromFile( char *filepath, double m[2000][2000] )
{
	FILE *inFile;
	ssize_t nread;
	size_t len = 0;
	char *line = NULL;
	char *token;
	int i = 0;
  int j = 0;
	inFile = fopen( filepath, "r" );

	while(( nread = getline( &line, &len, inFile )) != -1 ) {
		j = 0;
    token = strtok( line, " " );
		while( token != NULL && i < 2000 && j < 2000) {
		  m[i][j++] = atof( token );
			token = strtok( NULL, " " );
		}
    i++;
	}
	//free( line );
	//free( token );
	fclose( inFile );
}

int checkAnswer( int n, double result[2000][2000], double answer[2000][2000], int debug )
{
	int i, j;
	int correct = 1;

	for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
		  if ( fabs( result[i][j] - result[i][j] ) > EPSILON ) {
			  printf( "multiply answer: %.14lf	correct answer: %.14lf\n", result[i][j], answer[i][j]);
			  correct = 0;
		  } else if ( debug )  {
			  printf( "multiply answer: %.14lf	correct answer: %.14lf\n", result[i][j], answer[i][j]);
		  }
    }
	}
	return correct;
}
		
 		
void generateMatrix( int n, double *m )
{
	        int i;
		for ( i = 0; i < n * n; i++ ) {
			*(m + i) = rand() % 50; // mod 50 just for testing purposes
	        }

}
