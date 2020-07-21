#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h>
#include <inttypes.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/ioctl.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <pthread.h>

#include "matrixUtil.h"
#include "/home/valen/msr-safe/msr_safe.h"

#define N 2000
#define delta 0.05
#define NUM_THREADS_INDEX 1
#define PRINTED_FREQ 2.6 // in GHz

#define DEBUG 0
#define PRINT_MSRS 0
#define PRINT_FREQ 1



enum {
  BARRIER_INIT = 0,
  BARRIER_DELTA = 1,
  NUM_BARRIERS = 2
};

// struct to pass into threads
struct thread_info {
  int t;
  int start;
  int stop;
  double thread_start_time;
  double *calc_start_time;
  double *calc_stop_time;
  uint32_t grids;
  struct msr_batch_array *batch;
  struct timeval *msr_read_time;
  int num_threads;
};

// struct to hold clock frequency data
struct thread_msr_freq {
  double aperf_delta;
  double mperf_delta;
  double calc_freq;
  double mperf_time;
  double aperf_time;
};

int delta_result = 0;

int check_delta( double new[N][N], double old[N][N] ); 

void printGrid( double p[N][N] ); 

void initializeGrid( double p[N][N] );

void read_msrs( int threads, int fd, struct msr_batch_array *batch ); 

static pthread_barrier_t barrier[ NUM_BARRIERS ];

double a[N][N], b[N][N], correct[N][N];

// ************************** JACOBI ***************************
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
  if ( start_j == 0 || end_j == N ) {
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

// *************************** THREAD LOOP ***********************
void *thread_loop(void *threadnum) {
  struct thread_info *args  = ( struct thread_info *) (threadnum);
  int t = args->t;
  int start = args->start;
  int stop = args->stop;

  int count = 0;
  struct timeval init_start, init_stop, delta_start, delta_stop, calc_start, calc_stop, thread_start;
  double elapsed_delta= 0.0;
  double elapsed_calc = 0.0;
  
  //gettimeofday( &thread_start, NULL );
  //args->thread_start_time = thread_start.tv_sec + thread_start.tv_usec / 1000000.0;

  if (t == 0) {
    gettimeofday( &init_start, NULL);
	  initializeGrid( a );
	  initializeGrid( b );
    //gettimeofday( &init_stop, NULL);
    //fprintf(stdout, "%lf ", (init_stop.tv_sec - init_start.tv_sec) + (init_stop.tv_usec - init_start.tv_usec) / 1000000.0);
    //read_msrs( args->num_threads, args->batch, args->msr_read_time );
  }
  pthread_barrier_wait ( &barrier[BARRIER_INIT] );


  while(1) {
    //gettimeofday( &calc_start, NULL );
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
    //gettimeofday( &calc_stop, NULL );
    if (DEBUG) {
      args->calc_start_time[ count ] = ( calc_start.tv_sec - thread_start.tv_sec ) + ( calc_start.tv_usec - thread_start.tv_usec ) / 1000000.0;
      args->calc_stop_time[ count ] = ( calc_stop.tv_sec - thread_start.tv_sec ) + ( calc_stop.tv_usec - thread_start.tv_usec ) / 1000000.0;
    }
    count++;
    elapsed_calc += (calc_stop.tv_sec - calc_start.tv_sec) + (calc_stop.tv_usec - calc_start.tv_usec) / 1000000.0;

    //gettimeofday( &delta_start, NULL );
    if( t == 0) {
      delta_result = check_delta( a, b );
    }
    //gettimeofday( &delta_stop, NULL );
    elapsed_delta += (delta_stop.tv_sec - delta_start.tv_sec) + (delta_stop.tv_usec - delta_start.tv_usec) / 1000000.0;

    pthread_barrier_wait( &barrier[ BARRIER_DELTA ] );

    if( delta_result == 1 ) {
      break;
    }

  }

  if( t == 0 ) {
    ;
    //fprintf( stdout, "%lf %lf ", elapsed_delta, elapsed_calc );
  }

  //fprintf( stdout, "I am thread %d and elapsed calc is %lf\n", t,  elapsed_calc);
  pthread_exit( NULL );

}  

// *************************** MSR CODE ***************************	
void read_msrs( int threads, int fd, struct msr_batch_array *batch ) {
  int i, rc;

  for (i = 0; i < threads * 2; i+=2 ) {
    // mperf
    batch->ops[i].cpu = i/2;
    batch->ops[i].isrdmsr = 1;
    batch->ops[i].err = 0;
    batch->ops[i].msr = 0xE7;
    batch->ops[i].msrdata = 0;
    batch->ops[i].wmask = 0;

    // aperf
    batch->ops[i+1].cpu = i/2;
    batch->ops[i+1].isrdmsr = 1;
    batch->ops[i+1].err = 0;
    batch->ops[i+1].msr = 0xE8;
    batch->ops[i+1].msrdata = 0;
    batch->ops[i+1].wmask = 0;
  }


  //gettimeofday( read_time, NULL );
  rc = ioctl( fd, X86_IOC_MSR_BATCH, batch );
  assert( rc != -1 );
}  
		
void write_msrs( int threads, int fd, struct msr_batch_array *batch ) {
  int i, rc;

  for (i = 0; i < threads; i++ ) {
    // perf_ctl
    batch->ops[i].cpu = i/2;
    batch->ops[i].isrdmsr = 1;
    batch->ops[i].err = 0;
    batch->ops[i].msr = 0x199;
    batch->ops[i].msrdata = 0; // change this
    batch->ops[i].wmask = 0; // and this
  }


  //gettimeofday( read_time, NULL );
  rc = ioctl( fd, X86_IOC_MSR_BATCH, batch );
  assert( rc != -1 );
}  

void print_msrs( int threads, struct msr_batch_op start_op[], struct msr_batch_op stop_op[] ) {
  int i;
  for (i = 0; i < threads * 2; i+=2) {
    printf( "%2d mperf: %" PRIu64 "  %" PRIu64"  delta %" PRIu64"\n", i/2, (uint64_t)start_op[i].msrdata, (uint64_t)stop_op[i].msrdata, (uint64_t)stop_op[i].msrdata - (uint64_t)start_op[i].msrdata );
    printf( "   aperf: %" PRIu64 "  %" PRIu64"  delta %" PRIu64"\n", (uint64_t)start_op[i+1].msrdata, (uint64_t)stop_op[i+1].msrdata, (uint64_t)stop_op[i+1].msrdata - (uint64_t)start_op[i+1].msrdata );
  }
}

void calc_msr_freq( int threads, struct msr_batch_op start_op[], struct msr_batch_op stop_op[], struct thread_msr_freq answer[] ) {
  int i;
  for ( i = 0; i < threads * 2; i+=2 ) {
    answer[i/2].mperf_delta = ( (uint64_t) stop_op[i].msrdata - (uint64_t) start_op[i].msrdata );
    answer[i/2].aperf_delta = ( (uint64_t) stop_op[i+1].msrdata - (uint64_t) start_op[i+1].msrdata );
    answer[i/2].calc_freq = PRINTED_FREQ * ( answer[i/2].aperf_delta / answer[i/2].mperf_delta );
  }
}

//void print_thread_freq( int threads, struct thread_msr_freq freq[] ) {
//  int i;
//  for ( i = 0; i < threads; i++ ) {
//    printf( "thread %d   mperf freq: %.4lf  aperf freq: %.4lf\n", i, freq[i].mperf_freq, freq[i].aperf_freq );
//  }
//}


void print_thread_freq( int threads, struct thread_msr_freq freq[] ) {
  int i;
  for ( i = 0; i < threads; i++ ) {
    printf( "thread %2d   mperf delta: %.0lf  aperf delta: %.0lf  freq: %.4lf\n", i, freq[i].mperf_delta, freq[i].aperf_delta, freq[i].calc_freq);
  }
}

void calc_times( int threads, struct thread_msr_freq freq[] ) {
  int i;
  for( i = 0; i < threads; i++ ) {
    freq[i].mperf_time = freq[i].mperf_delta / PRINTED_FREQ / 1.0E9;
    freq[i].aperf_time = freq[i].aperf_delta / freq[i].calc_freq / 1.0E9;
  }
}

void print_thread_times( int threads, struct thread_msr_freq freq[] ) {
  int i;
  for( i = 0; i < threads; i++ ) {
    printf( "thread %2d   mperf time: %.6lf  aperf time: %.6lf\n", i, freq[i].mperf_time, freq[i].aperf_time );
  }
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
  int num_threads, t, step, fd;

  double elapsed_time, begin_time, end_time;
  
  num_threads = atoi( argv[ NUM_THREADS_INDEX ] );

  struct msr_batch_array batch;
  struct msr_batch_op start_op[ num_threads * 2 ], after_init_op[ num_threads * 2], stop_op[ num_threads * 2];
  struct thread_msr_freq freq_total[ num_threads ], freq_begin[ num_threads ], freq_end[ num_threads ];

  struct timeval msr_read_start, after_init_time, msr_read_stop;
  batch.numops = num_threads * 2;
  batch.ops = start_op;

  fd = open( "/dev/cpu/msr_batch", O_RDWR );

  read_msrs( num_threads, fd, &batch );

  step = ceil( N / ( double ) num_threads );

  struct thread_info thread_args[ num_threads ];
  pthread_t threads[ num_threads ];

  for( t = 0; t < NUM_BARRIERS; t++ ) {
    assert( ! pthread_barrier_init( &barrier[t], NULL, num_threads ) );
  }

  //batch.ops = after_init_op;
  for( t = 0; t < num_threads; t++ ) {
    thread_args[ t ].t = t;
    thread_args[ t ].start = t * step;
    thread_args[ t ].stop = ( t != num_threads - 1 ) ? ( step * ( t + 1 ) ) : N;
    //thread_args[ t ].thread_start_time = 0;
    //thread_args[ t ].calc_start_time = ( double * ) malloc( sizeof( double ) * 200 );
    //thread_args[ t ].calc_stop_time = ( double * ) malloc( sizeof( double ) * 200 );
    //thread_args[ t ].batch = &batch;
    //thread_args[ t ].msr_read_time = &after_init_time;
    //thread_args[ t ].num_threads = num_threads;
    assert( ! pthread_create( &threads[ t ], NULL, thread_loop, (void*)&thread_args[t] ) );
  }

  for( t = 0; t < num_threads; t++ ) {
    pthread_join( threads[ t ], NULL );
  }

  batch.ops = stop_op;
  read_msrs( num_threads, fd, &batch );

  for( t = 0; t < NUM_BARRIERS; t++ ) {
    assert( ! pthread_barrier_destroy( &barrier[ t ] ) );
  }

  //elapsed_time = ( msr_read_stop.tv_sec - msr_read_start.tv_sec ) + ( msr_read_stop.tv_usec - msr_read_start.tv_usec) / 1000000.0;
  //begin_time = ( after_init_time.tv_sec - msr_read_start.tv_sec ) + ( after_init_time.tv_usec - msr_read_start.tv_usec) / 1000000.0;
  //end_time = ( msr_read_stop.tv_sec - after_init_time.tv_sec ) + ( msr_read_stop.tv_usec - after_init_time.tv_usec) / 1000000.0;

  //calc_msr_freq( num_threads, elapsed_time, start_op, stop_op, freq_total );
  //calc_msr_freq( num_threads, begin_time, start_op, after_init_op, freq_begin );
  //calc_msr_freq( num_threads, end_time, after_init_op, stop_op, freq_end );
  calc_msr_freq( num_threads, start_op, stop_op, freq_total );

  if ( PRINT_FREQ ) {
    printf( "\n" );
    //print_thread_freq( num_threads, freq );
    print_thread_freq( num_threads, freq_total );
    printf( "\n" );
    calc_times( num_threads, freq_total );
    print_thread_times( num_threads, freq_total );
  }

  if (DEBUG) {
    int i;
    printf("\n");
    for ( t = 0; t < num_threads; t++ ) {
      printf("%10.6f ", thread_args[ t ].thread_start_time);
    }
    printf("\n");
    for ( i = 0; i < 200; i++ ) {
      for ( t = 0; t < num_threads; t++ ) {
        printf("%10.6lf ", thread_args[ t ].calc_start_time[ i ] );
      }
      printf("\n");
      for ( t = 0; t < num_threads; t++ ) {
        printf("%10.6lf ", thread_args[ t ].calc_stop_time[ i ] );
      }
      printf("\n\n\n");

    }

  }

  if (PRINT_MSRS) {
    printf("\n");
    print_msrs( num_threads, start_op, stop_op );
    //printf( "Elapsed time: %.4lf\n", elapsed_time );
    //printf( "Begin time: %.4lf\n", begin_time );
    //printf( "End time: %.4lf\n", end_time );
  }
  //readMatrixFromFile( "correct", correct );
  //checkAnswer( N, b, correct, 0);

  close( fd );
  pthread_exit(NULL);

}
	
