#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <pthread.h>

#include "matrixUtil.h"

#define N 2000
#define delta 0.05
#define NUM_THREADS_INDEX 1


enum {
  BARRIER_INIT = 0,
  BARRIER_DELTA = 1,
  NUM_BARRIERS = 2
};

struct thread_info {
  int t;
  int start;
  int stop;
};

int delta_result = 0;

int check_delta( double new[N][N], double old[N][N] ); 

void printGrid( double p[N][N] ); 

void initializeGrid( double p[N][N] );

static pthread_barrier_t barrier[ NUM_BARRIERS ];

double a[N][N], b[N][N], correct[N][N];

void jacobi( double p[N][N], double q[N][N], int start_j, int end_j ) {
	int i, j;
  int end = ( end_j == N ) ? N - 1 : end_j;
  // middle 
	for ( i = 1; i < N - 1; i++ ) {
    j = (start_j == 0) ? 1: start_j;
		for ( ; j < end; j++) {
			q[i][j] = p[i+1][j-1] + p[i+1][j] + p[i+1][j+1]
				+ p[i][j-1] + p[i][j] + p[i][j+1]
				+ p[i-1][j-1] + p[i-1][j] + p[i-1][j+1];
			q[i][j] /= 9.0;
		}
	}
  // vertical sides
	for ( i = 1; i < N-1; i ++) {
    if ( start_j == 0 ){
		  q[i][0] = p[i+1][0] + p[i+1][1]
			  + p[i][0] + p[i][1]	
			  + p[i-1][0] + p[i-1][1];
		  q[i][0] /= 6.0;
    }
    if ( end_j == N ) {
		  q[i][N-1] = p[i+1][N-2] + p[i+1][N-1] 
			  + p[i][N-2] + p[i][N-1] 
			  + p[i-1][N-2] + p[i-1][N-1];
		  q[i][N-1] /=6.0;
    }
	}
  // horizonal sides
  j = (start_j == 0) ? 1: start_j;
	for ( ; j < end; j++) {
		i = 0;
		q[0][j] = p[i][j-1] + p[i][j] + p[i][j+1]
			+ p[i+1][j-1] + p[i+1][j] + p[i+1][j+1];

		q[0][j] /= 6.0;
		i = N- 1;
		q[N-1][j] = p[i-1][j-1] + p[i-1][j] + p[i-1][j+1]
			+ p[i][j-1] + p[i][j] + p[i][j+1];
		q[N-1][j] /= 6.0;
	}
  // two corners
  if (start_j == 0 ) {
	  q[0][N-1] = p[0][N-1] + p[0][N-2] + p[1][N-2] + p[1][N-1];
	  q[0][N-1] /= 4.0;
  }
  if ( end_j == N ) {
	  q[N-1][0] = p[N-1][0] + p[N-2][0] + p[N-2][1] + p[N-1][1];
	  q[N-1][0] /= 4.0;
  }
} 

void *thread_loop(void *threadnum) {
  struct thread_info *args  = ( struct thread_info *) (threadnum);
  int t = args->t;
  int start = args->start;
  int stop = args->stop;

  int count = 0;
  struct timeval init_start, init_stop, delta_start, delta_stop, calc_start, calc_stop;
  double elapsed_delta= 0.0;
  double elapsed_calc = 0.0;

  if (t == 0) {
    gettimeofday( &init_start, NULL);
	  initializeGrid( a );
	  initializeGrid( b );
    gettimeofday( &init_stop, NULL);
    fprintf(stdout, "%lf ", (init_stop.tv_sec - init_start.tv_sec) + (init_stop.tv_usec - init_start.tv_usec) / 1000000.0);
  }
  pthread_barrier_wait ( &barrier[BARRIER_INIT] );

  while(1) {
    count++;
    gettimeofday( &calc_start, NULL );
    if (! (count%2) ) {
      jacobi( a, b, start, stop);
//      if ( t == 0) {
//        printf( " B\n" );
//        printGrid(b);
//      }
    } else {
      jacobi( b, a, start, stop);
//      if( t == 0 ) {
//        printf( " A\n" );
//        printGrid(a);
//      }
    }
    gettimeofday( &calc_stop, NULL );
    elapsed_calc += (calc_stop.tv_sec - calc_start.tv_sec) + (calc_stop.tv_usec - calc_start.tv_usec) / 1000000.0;

    gettimeofday( &delta_start, NULL );
    if( t == 0) {
      delta_result = check_delta( a, b );
    }
    gettimeofday( &delta_stop, NULL );
    elapsed_delta += (delta_stop.tv_sec - delta_start.tv_sec) + (delta_stop.tv_usec - delta_start.tv_usec) / 1000000.0;

    pthread_barrier_wait( &barrier[ BARRIER_DELTA ] );

    if( delta_result == 1 ) {
      break;
    }

  }

  if( t == 0 ) {
    fprintf( stdout, "%lf %lf\n", elapsed_delta, elapsed_calc );
  }
  pthread_exit( NULL );

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

int main( int argc, char *argv[] )  {
  int num_threads, t, step;
  
  num_threads = atoi( argv[ NUM_THREADS_INDEX ] );

  step = ceil( N / num_threads );

  struct thread_info thread_args[ num_threads ];
  pthread_t threads[ num_threads ];

  for( t = 0; t < NUM_BARRIERS; t++ ) {
    assert( ! pthread_barrier_init( &barrier[t], NULL, num_threads ) );
  }

  for( t = 0; t < num_threads; t++ ) {
    thread_args[ t ].t = t;
    thread_args[ t ].start = t * step;
    thread_args[ t ].stop = ( t != num_threads - 1 ) ? ( step * ( t + 1 ) ) : N;
    assert( ! pthread_create( &threads[ t ], NULL, thread_loop, (void*)&thread_args[t] ) );
  }

  for( t = 0; t < num_threads; t++ ) {
    pthread_join( threads[ t ], NULL );
  }

  for( t = 0; t < NUM_BARRIERS; t++ ) {
    assert( ! pthread_barrier_destroy( &barrier[ t ] ) );
  }

  readMatrixFromFile( "correct", correct );
  checkAnswer( N, b, correct, 0);
  pthread_exit(NULL);

}
	
