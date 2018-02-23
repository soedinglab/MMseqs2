#!/bin/bash -e
REPO="$(readlink -f $1)"
BUILD="$(readlink -f $2)"

if [ ! -d "$REPO" ]; then
    echo "MMseqs2 repository missing"
    exit 1
fi

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
mkdir -p "$BUILD/mmseqs/bin"
objcopy src/mmseqs --compress-debug-sections
strip --only-keep-debug src/mmseqs -o "$BUILD/mmseqs_debug_symbols_sse41"
strip --strip-debug src/mmseqs
cp src/mmseqs "$BUILD/mmseqs/bin/mmseqs_sse41"
for i in $(ldd src/mmseqs | awk '{ print $3 }' | grep -v cygdrive | grep -v '???'); do
    cp $i "$BUILD/mmseqs/bin";
done

cd "$BUILD" && mkdir -p "build_avx2" && cd "build_avx2"
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_AVX2=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
objcopy src/mmseqs --compress-debug-sections
strip --only-keep-debug src/mmseqs -o "$BUILD/mmseqs_debug_symbols_avx2"
strip --strip-debug src/mmseqs
cp src/mmseqs "$BUILD/mmseqs/bin/mmseqs_avx2"

cp -f /usr/libexec/busybox-standalone/bin/busybox.exe "$BUILD/mmseqs/bin"
cd "$BUILD/mmseqs/bin"
cmd /c ""busybox.exe" "--install" "-s" ".""

cat <<'CPUDOC' | gcc -Os -std=gnu99 -o "$BUILD/mmseqs/bin/testcpu.exe" -xc -
#include <stdio.h>
#define D(x) __builtin_cpu_supports((x))?(x)
int main(){puts(D("avx2"):D("sse4.1"):"fail");}
CPUDOC

cat <<'BATDOC' > "$BUILD/mmseqs/mmseqs.bat"
for /f %%i in ('%~dp0\bin\testcpu.exe') do (
if "%%i" == "avx2" ( %~dp0\bin\mmseqs_avx2.exe %* )
if "%%i" == "sse4.1" ( %~dp0\bin\mmseqs_sse41.exe %* )
if "%%i" == "fail" ( echo Unsupported CPU! MMseqs2 requires SSE4.1 or AVX2 instruction set support. )
)
BATDOC

chmod +x "$BUILD/mmseqs/mmseqs.bat"
