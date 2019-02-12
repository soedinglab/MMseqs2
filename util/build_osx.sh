#!/bin/bash -e
REPO="$(greadlink -f $1)"
BUILD="$(greadlink -f $2)"

if [ ! -d "$REPO" ]; then
    echo "MMseqs2 repository missing"
    exit 1
fi

mkdir -p "$BUILD/mmseqs"
mkdir -p "$BUILD/lib"
cd "$BUILD/lib"
ln -s `$CXX --print-file-name=libgomp.a`
export CXXFLAGS="-L$BUILD/lib"

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$BUILD/mmseqs" -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
make install

cd "$BUILD"
tar -czvf mmseqs-osx-static_sse41.tar.gz mmseqs

cd "$BUILD" && mkdir -p "build_avx2" && cd "build_avx2"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_AVX2=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
cp src/mmseqs "$BUILD/mmseqs/bin/mmseqs"

cd "$BUILD"
tar -czvf mmseqs-osx-static_avx2.tar.gz mmseqs
