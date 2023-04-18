#!/usr/bin/bash

DIR=$1

if [ -z ${DIR} ]
then
  echo "Usage: $0 <model directory>"
  exit 0
fi

# Cleaning

find . -name "*.nc" -delete
find . -name "*.log" -delete
find . -name "SFINCStest-*" -delete

# Running qsub tasks

for MODEL in `\ls $DIR`
do
  echo "Running model in $MODEL"
  qsub -N SFINCStest-${MODEL} test-job.sh ${DIR} ${MODEL}
done

exit

