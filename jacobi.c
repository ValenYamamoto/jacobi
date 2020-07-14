#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#define N 20
#define delta 0.1


void jacobi( double p[N][N], double q[N][N]) {
	int i, j;
	for ( i = 1; i < N - 1; i++ ) {
		for ( j = 1; j < N-1; j++) {
			q[i][j] = p[i+1][j-1] + p[i+1][j] + p[i+1][j+1]
				+ p[i][j-1] + p[i][j] + p[i][j+1]
				+ p[i-1][j-1] + p[i-1][j] + p[i-1][j+1];
			q[i][j] /= 9.0;
		}
	}
	for ( i = 1; i < N-1; i ++) {
		q[i][0] = p[i+1][0] + p[i+1][1]
			+ p[i][0] + p[i][1]	
			+ p[i-1][0] + p[i-1][1];
		q[i][0] /= 6.0;
		q[i][N-1] = p[i+1][N-2] + p[i+1][N-1] 
			+ p[i][N-2] + p[i][N-1] 
			+ p[i-1][N-2] + p[i-1][N-1];
		q[i][N-1] /=6.0;
	}
	for ( j = 1; j < N-1; j++) {
		i = 0;
		q[0][j] = p[i][j-1] + p[i][j] + p[i][j+1]
			+ p[i+1][j-1] + p[i+1][j] + p[i+1][j+1];

		q[0][j] /= 6.0;
		i = N- 1;
		q[N-1][j] = p[i-1][j-1] + p[i-1][j] + p[i-1][j+1]
			+ p[i][j-1] + p[i][j] + p[i][j+1];
		q[N-1][j] /= 6.0;
	}
	q[0][N-1] = p[0][N-1] + p[0][N-2] + p[1][N-2] + p[1][N-1];
	q[0][N-1] /= 4.0;
	q[N-1][0] = p[N-1][0] + p[N-2][0] + p[N-2][1] + p[N-1][1];
	q[N-1][0] /= 4.0;
} 

	

	
	
		

void initializeGrid( double p[N][N] ) {
	int i, j;
	for ( i = 0; i < N; i ++ ) {
		for ( j= 0; j < N; j++) {
			p[i][j] = 0.0;
		}
	}
	p[0][0] = -100.0;
	p[N-1][N-1] = 100.0;
}

void printGrid( double p[N][N] ) {
	int i, j;
	printf("\n");
	for ( i = 0; i < N; i ++ ) {
		for ( j= 0; j < N; j++ ) {
			printf( "%6.2f  ", p[i][j]);
		}
		printf("\n");
	}
}

int check_delta( double new[N][N], double old[N][N] ) {
	int i, j; 
	double max_delta = 0;
	for ( i = 0; i < N; i ++ ) {
		for ( j = 0; j < N; j ++ ) {
			if  ( fabs( new[i][j] - old[i][j] ) > max_delta ) {
				max_delta =  fabs( new[i][j] - old[i][j] );
			}
		}
	}
	if ( max_delta < delta) {
		return 1;
	}
	return 0;
}

int main() {
	double a[N][N], b[N][N];
	int exetime;


	struct timeval start, end;

	initializeGrid( a );
	initializeGrid( b );
	
	gettimeofday( &start, NULL );	
	while (1) {
		jacobi(a, b);
//		printGrid( b);
		if ( check_delta(b, a) == 1) {
			gettimeofday( &end, NULL );
//			printGrid(b);
			break;
		}
		jacobi(b, a);
//		printGrid(a);
		if ( check_delta(a, b) == 1 ) {
			gettimeofday( &end, NULL );
//			printGrid(a);
			break;
		}
	}	
	exetime = ( end.tv_sec * 1000000 + end.tv_usec ) - ( start.tv_sec * 1000000 + start.tv_usec );
	printf( "%d", exetime );
	return 0;
}
	
