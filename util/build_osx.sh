#!/bin/sh -e
abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

REPO="$(abspath "$1")"
BUILD="$(abspath "$2")"
BINARY_NAME="${3:-mmseqs}"

if [ ! -d "$REPO" ]; then
    echo "${BINARY_NAME} repository missing"
    exit 1
fi

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/usr/local/opt/libomp/lib/libomp.a "$REPO"
make -j 4

if [ "$(echo $(otool -L "src/${BINARY_NAME}" | wc -l))" != 5 ]; then
    echo "Too many linked libraries found in ${BINARY_NAME} binary. Build is not static!"
    exit 1
fi

cd "$BUILD" && mkdir -p "build_avx2" && cd "build_avx2"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_AVX2=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/usr/local/opt/libomp/lib/libomp.a "$REPO"
make -j 4

if [ "$(echo $(otool -L "src/${BINARY_NAME}" | wc -l))" != 5 ]; then
    echo "Too many linked libraries found in ${BINARY_NAME} binary. Build i not static!"
    exit 1
fi

