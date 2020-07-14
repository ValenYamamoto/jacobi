#define EPSILON 1E-7


void readMatrixFromFile( char *filepath, double m[2000][2000] );

int checkAnswer( int n, double result[2000][2000], double answer[2000][2000], int debug );

void generateMatrix( int n, double *m );
