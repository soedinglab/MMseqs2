#!/bin/bash -x
docker build -f Dockerfile.build -t mmseqs-build-image .
docker run -t mmseqs-build-image:latest
docker cp $(docker ps -l -q):/root/mmseqs/ .
docker build -f Dockerfile.release -t mmseqs .
CURR_BUILD=$(pwd)/mmseqs-$(date +%Y-%m-%d-%T)
mkdir -p $CURR_BUILD
cp -r mmseqs/mmseqs/bin $CURR_BUILD
cp -r mmseqs/data $CURR_BUILD
cp -r mmseqs/util $CURR_BUILD
rm -rf $CURR_BUILD/util/docker
cp -r mmseqs/userguide.pdf $CURR_BUILD
tar cvfz $(pwd)/mmseqs-static.tar.gz -C $CURR_BUILD bin data util userguide.pdf
mv mmseqs $CURR_BUILD
docker run mmseqs 
