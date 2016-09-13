#!/bin/bash -e
cd /root/mmseqs
git checkout .
git clean -f -d
git pull
cd /root/mmseqs/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/root/mmseqs/mmseqs/ -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS_RELEASE="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" /root/mmseqs
make VERBOSE=1
make install
