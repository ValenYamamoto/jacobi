#!/bin/python3
import os
import matplotlib.pyplot as plt
from collections import defaultdict
from math import ceil

def process_line( line: str ):
  processed = tuple(x  for x in line.rstrip( "\n" ).split( " " ) ) 
  init_time, delta_time, calc_time, real_time, user_time = tuple( float( x ) for x in processed[:-1] )
  return init_time, delta_time, calc_time, real_time, user_time

if __name__ == "__main__":
  times = defaultdict( list )

  with open( "./test" ) as f:
    for num, line in enumerate( f, 1 ):
      print(line);
      init_time, delta_time, calc_time, real_time, user_time = process_line( line )
      times[ ceil( num / 3 ) ] += [ real_time ]

  for threads, exectimes in times.items():
    plt.scatter( [threads] * 3, exectimes )
  plt.savefig( "graph_real.png" )
