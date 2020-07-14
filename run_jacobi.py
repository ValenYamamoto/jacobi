import os
import sys
import numpy

def process_stream( stream: "iterable" ) -> int:
	iterable = iter( stream )
	return float( next( iterable ).rstrip( "\n" ) )

if __name__ == "__main__":
	answers = []
	for _ in range( 100 ):
		stream = os.popen( './a.out' )
		exetime = process_stream( stream )
		answers.append( exetime )

	print( f"mean: {numpy.mean( answers )} " )
	print( f"std dev: {numpy.std( answers, ddof=1 )} " )
