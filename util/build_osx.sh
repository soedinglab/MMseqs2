#!/bin/bash -e
REPO="$(greadlink -f $1)"
BUILD="$(greadlink -f $2)"
BINARY_NAME="${3:-mmseqs}"

if [ ! -d "$REPO" ]; then
    echo "${BINARY_NAME} repository missing"
    exit 1
fi

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4

if [ "$(echo $(otool -L "src/${BINARY_NAME}" | wc -l))" != 4 ]; then
    echo "Too many linked libraries found in ${BINARY_NAME} binary. Build is not static!"
    exit 1
fi

cd "$BUILD" && mkdir -p "build_avx2" && cd "build_avx2"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_AVX2=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4

if [ "$(echo $(otool -L "src/${BINARY_NAME}" | wc -l))" != 4 ]; then
    echo "Too many linked libraries found in ${BINARY_NAME} binary. Build i not static!"
    exit 1
fi

