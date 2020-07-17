#!/bin/bash

EXEC_FILE=$1
NUM_THREADS=$2
for (( THREADS=1; THREADS<=NUM_THREADS; THREADS++ )); do
  for (( i=0; i<3; i++ )); do
    /usr/bin/time --format "%e %U" "./$EXEC_FILE" $THREADS >> test 2>&1
  done
done
