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

export MACOSX_DEPLOYMENT_TARGET=10.11

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/usr/local/opt/libomp/lib/libomp.a "$REPO"
make -j $(nproc --all)

if [ "$(echo $(otool -L "src/${BINARY_NAME}" | wc -l))" != 5 ]; then
    echo "Too many linked libraries found in ${BINARY_NAME} binary. Build is not static!"
    exit 1
fi

cd "$BUILD" && mkdir -p "build_avx2" && cd "build_avx2"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_AVX2=1 -DCMAKE_CXX_FLAGS="-arch x86_64h" -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/usr/local/opt/libomp/lib/libomp.a "$REPO"
make -j $(nproc --all)

if [ "$(echo $(otool -L "src/${BINARY_NAME}" | wc -l))" != 5 ]; then
    echo "Too many linked libraries found in ${BINARY_NAME} binary. Build is not static!"
    exit 1
fi

if sw_vers -productVersion | awk -F. '$1 < 11 && $2 < 16 { exit 0; } { exit 1; }'; then
    lipo -create -arch x86_64 "$BUILD/build_sse41/src/${BINARY_NAME}" -arch x86_64h "$BUILD/build_avx2/src/${BINARY_NAME}" -output "$BUILD/${BINARY_NAME}"
    exit 0
fi

unset MACOSX_DEPLOYMENT_TARGET

cd "$BUILD" && mkdir -p "build_libomp" && cd "build_libomp"
wget -qO- http://github.com/llvm/llvm-project/releases/download/llvmorg-11.0.0/openmp-11.0.0.src.tar.xz | tar xvf -
cd openmp-11.0.0.src
wget https://raw.githubusercontent.com/Homebrew/formula-patches/7e2ee1d7/libomp/arm.patch
patch -p1 < arm.patch
mkdir build && cd build
cmake -DLIBOMP_ENABLE_SHARED=OFF -DLIBOMP_INSTALL_ALIASES=OFF -DLIBOMP_ARCH=aarch64 -DCMAKE_CXX_FLAGS="-arch arm64 -target arm64-apple-macos11" -DLIBOMP_ASMFLAGS="-arch arm64 -target arm64-apple-macos11" ..
make -j $(nproc --all)
export LIBOMP_AARCH64="$BUILD/build_libomp/openmp-11.0.0.src/build/runtime/src"

cd "$BUILD" && mkdir -p "build_arm64" && cd "build_arm64"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_ARM8=1 -DCMAKE_C_FLAGS="-arch arm64 -target arm64-apple-macos11" -DCMAKE_CXX_FLAGS="-arch arm64 -target arm64-apple-macos11" -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AARCH64}" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AARCH64}" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=${LIBOMP_AARCH64}/libomp.a "$REPO"
make -j $(nproc --all)
if [ "$(echo $(otool -L "src/${BINARY_NAME}" | wc -l))" != 5 ]; then
    echo "Too many linked libraries found in ${BINARY_NAME} binary. Build is not static!"
    exit 1
fi

lipo -create -arch x86_64 "$BUILD/build_sse41/src/${BINARY_NAME}" -arch x86_64h "$BUILD/build_avx2/src/${BINARY_NAME}" -arch arm64 "$BUILD/build_arm64/src/${BINARY_NAME}" -output "$BUILD/${BINARY_NAME}"
