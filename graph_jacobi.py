#!/bin/python3
import os
import matplotlib.pyplot as plt
from collections import defaultdict

def process_stream( stream: "iterable" ) -> int:
  output = stream.read()
  print(output)
  init_time, delta_time, calc_time = tuple( float( x ) for x in output.rstrip( "\n" ).split( " " ) )
  return init_time, delta_time, calc_time

if __name__ == "__main__":
  times = defaultdict( list )

  for i in range( 1, 100 ): # number of cores
    for _ in range(3):
      stream = os.popen( f'./a.out {i}' )
      init_time, delta_time, calc_time = process_stream( stream )
      stream.close()
      times[i].append( calc_time )

  for threads, exectimes in times.items():
    plt.scatter( [threads] * 3, exectimes )
  plt.savefig( "graph100.png" )
