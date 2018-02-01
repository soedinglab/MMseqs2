#!/bin/bash -e
REPO="$(greadlink -f $1)"
BUILD="$(greadlink -f $2)"

if [ ! -d "$REPO" ]; then
    echo "MMseqs2 repository missing"
    exit 1
fi

mkdir -p "$BUILD/mmseqs"

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
# TODO FIXME: Static build not possible with gcc on MacOS?
# complains about missing lcrt, remove _INVALID suffix here and further down once we find a solution
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="$BUILD/mmseqs" -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS_INVALID="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
make install

gobjcopy "$BUILD/mmseqs/bin/mmseqs" --compress-debug-sections
gstrip --only-keep-debug "$BUILD/mmseqs/bin/mmseqs" -o "$BUILD/mmseqs_debug_symbols_sse41"
gstrip --strip-debug "$BUILD/mmseqs/bin/mmseqs"
cd "$BUILD"
tar -czvf mmseqs-osx-static_sse41.tar.gz mmseqs

cd "$BUILD" && mkdir -p "build_avx2" && cd "build_avx2"
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_AVX2=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS_INVALID="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
gobjcopy src/mmseqs --compress-debug-sections
gstrip --only-keep-debug src/mmseqs -o "$BUILD/mmseqs_debug_symbols_avx2"
gstrip --strip-debug src/mmseqs
cp src/mmseqs "$BUILD/mmseqs/bin/mmseqs"

cd "$BUILD"
tar -czvf mmseqs-osx-static_avx2.tar.gz mmseqs
tar -czvf mmseqs-osx-debug-symbols.tar.gz mmseqs_debug_symbols_sse41 mmseqs_debug_symbols_avx2
