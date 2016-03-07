#!/bin/bash -x
docker run -t mmseqs-build-image:latest
docker cp $(docker ps -l -q):/root/mmseqs/ .
docker build -f Dockerfile.release -t mmseqs .
CURR_BUILD=$(pwd)/mmseqs-$(date +%Y-%m-%d-%T)
mkdir -p $CURR_BUILD
cp -r mmseqs/bin/ $CURR_BUILD
cp -r mmseqs/data/ $CURR_BUILD
cp -r mmseqs/userguide.pdf $CURR_BUILD
tar cvfz $(pwd)/mmseqs-static.tar.gz -C $CURR_BUILD bin data userguide.pdf
mv mmdir $CURR_BUILD
docker run mmseqs 
