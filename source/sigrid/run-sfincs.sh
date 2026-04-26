#!/usr/bin/bash

IMPORT_DIRECTORY=$1

if [ -z ${IMPORT_DIRECTORY} ]
then
  echo "Usage: $0 <directory in container which contains input>"
  echo "For example: $0 /data if the Singularity container is started as: ./sfincs-cpu --bind <input>:/data"
  exit 0
fi

cd ${IMPORT_DIRECTORY}
/usr/local/bin/sfincs

exit
