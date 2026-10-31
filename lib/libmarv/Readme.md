




# LIBMARV

GPU-accelerated protein database scan.

## Software requirements
* (NVIDIA) Linux operating system with compatible CUDA Toolkit 13.2 or newer
* (AMD) Linux operating system with compatible ROCm 7.2.4 or newer
* C++17 compiler
* zlib
* make
* cmake 3.24

## Hardware requirements
* (NVIDIA) A modern CUDA-capable GPU of Turing generation or newer. We provide tunings up to the Blackwell architecture. By PTX forward compatibility, libmarv should run on newer GPUs although performance will be less than expected without dedicated tuning.
* (AMD) A modern ROCm-capable GPU. Currently, we only provide tunings for gfx900, gfx942, gfx1100, and gfx1150, and do not offer a compatibility fallback for other architectures.




## Build

Libmarv is built using CMake.

```
mkdir build
cd build
cmake ../src/CMakeLists.txt -DCMAKE_BUILD_TYPE=Release -DENABLE_HIP=0 -DBUILD_TUNING=0 -DONLY_BUILD_TUNING=0 -DLIBRARY_ONLY=0 -DCMAKE_CUDA_ARCHITECTURES="120"
cmake ../src/CMakeLists.txt -DCMAKE_BUILD_TYPE=Release -DENABLE_HIP=0 -DBUILD_TUNING=0 -DONLY_BUILD_TUNING=0 -DCMAKE_CUDA_ARCHITECTURES="120"
```

Use `ENABLE_HIP=0` to target NVIDIA GPUs and `ENABLE_HIP=1` to target AMD GPUs. The cmake variables CMAKE_CUDA_ARCHITECTURES and CMAKE_HIP_ARCHITECTURES specify the GPU architectures to compile for NVIDIA or AMD, respectively. All architectures on which libmarv is intended to be used should be specified for best performance.

With `LIBRARY_ONLY=1` only the library part used by MMseqs2 will be built. `LIBRARY_ONLY=0` will additionally build a standalone program which can serve as an example on how to use the different database searches.

When `BUILD_TUNING=1` our tuning program for alignment kernels will be built together with the library part and standalone program. 
To only build the tuning program, pass `ONLY_BUILD_TUNING=1`

Example for building library and standalone program targetting NVIDIA H100 and NVIDIA RTX Pro 6000 Blackwell
```
cmake ../src/CMakeLists.txt -DCMAKE_BUILD_TYPE=Release -DENABLE_HIP=0 -DBUILD_TUNING=0 -DONLY_BUILD_TUNING=0 -DLIBRARY_ONLY=0 -DCMAKE_CUDA_ARCHITECTURES="90;120"
```


## Standalone program

The standalone program has two components, **makedb** and **align** . **makedb** is used to construct a database which can be queried by **align**.

## Database construction
Use **makedb** to create a database from a fasta file. The file can be gzip'ed.
We support fasta files with up to 2 billion sequences.

```
mkdir -p dbfolder
./makedb input.fa(.gz) dbfolder/dbname [options]
```

Options:
* --mem val : Memory limit. Can use suffix K,M,G. If makedb requires more memory, temp files in temp directory will be used. Default all available memory.
* --tempdir val : Temp directory for temporary files. Must exist. Default is db output directory.
* --allowLowerCase : Convert lower-case letters to distinct numbers. Without this option, lower-case letters are treated as upper-case.\n";



## Querying the database
Use **align** to query the database. The environment variables CUDA_VISIBLE_DEVICES / HIP_VISIBLE_DEVICES control which NVIDIA / AMD GPUs in the system should be used. If the corresponding environment variable is not set, all GPUs will be used.

**align** has two mandatory arguments. 
1. `--query` The query file which contains all queries
2. `--db` The path to the reference database constructed with makedb. 

Run `./align --help` to get a complete list of options.

Example call:
```
./align --query queries.fa(.gz) --db dbfolder/dbname -prefetchDBFile --scanType Gapless --top 10 --maxGpuMem 27G --uploadFull --queryBatchsize 5
```


## Metal and SIMD backends

`-DENABLE_METAL=ON` builds the Apple GPU backend (macOS, Xcode Metal compiler), `-DENABLE_SIMD=ON` the CPU backend
(NEON, SSE2, AVX2, AVX-512BW). Both need a C++20 compiler and OpenMP, e.g.
`cmake -S src -B build -DENABLE_SIMD=ON -DCMAKE_BUILD_TYPE=Release -DOpenMP_ROOT=/opt/homebrew/opt/libomp`.
`Marv::MemoryConfig` bounds the database buffers per device (`maxGpuMem`) and the score-only workspace kept between
searches (`maxCachedQueryBytes`, 0 disables reuse); `releaseSearchWorkspace()` frees it explicitly.
