#!/bin/bash -e
docker build -f Dockerfile.build -t mmseqs-build-image .
docker run -t mmseqs-build-image:latest
docker cp $(docker ps -l -q):/root/mmseqs .
docker build -f Dockerfile.release -t mmseqs .
docker run mmseqs 
rm -rf mmseqs
