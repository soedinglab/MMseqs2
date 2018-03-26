#!/bin/bash -e
REPO="$(readlink -f $1)"
BUILD="$(readlink -f $2)"

if [ ! -d "$REPO" ]; then
    echo "MMseqs2 repository missing"
    exit 1
fi

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
mkdir -p "$BUILD/mmseqs/bin"
cp src/mmseqs "$BUILD/mmseqs/bin/mmseqs_sse41"
for i in $(ldd src/mmseqs | awk '{ print $3 }' | grep -v cygdrive | grep -v '???'); do
    cp $i "$BUILD/mmseqs/bin";
done

cd "$BUILD" && mkdir -p "build_avx2" && cd "build_avx2"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_AVX2=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
cp src/mmseqs "$BUILD/mmseqs/bin/mmseqs_avx2"

cp -f /usr/libexec/busybox-standalone/bin/busybox.exe "$BUILD/mmseqs/bin"

cat <<'CPUDOC' | gcc -Os -std=gnu99 -o "$BUILD/mmseqs/bin/testcpu.exe" -xc -
#include <stdio.h>
#define D(x) __builtin_cpu_supports((x))?(x)
int main(){puts(D("avx2"):D("sse4.1"):"fail");}
CPUDOC

cp "$REPO/util/mmseqs_wrapper.bat" "$BUILD/mmseqs/mmseqs.bat"
chmod +x "$BUILD/mmseqs/mmseqs.bat"
