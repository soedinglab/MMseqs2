ARG APP=mmseqs
FROM --platform=$BUILDPLATFORM debian:bookworm-slim AS builder
ARG TARGETARCH
ARG APP
ARG GPU

RUN dpkg --add-architecture $TARGETARCH \
    && apt-get update \
    && apt-get install -y \
      build-essential cmake xxd git wget \
      zlib1g-dev libbz2-dev libatomic1 \
      crossbuild-essential-$TARGETARCH zlib1g-dev:$TARGETARCH libbz2-dev:$TARGETARCH; \
    if [ "$GPU" = "1" ]; then \
      wget https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64/cuda-keyring_1.1-1_all.deb; \
      dpkg -i cuda-keyring_1.1-1_all.deb; \
      apt-get update && apt-get install -y cuda-nvcc-12-6 cuda-cudart-dev-12-6 ninja-build; \
    fi; \
    rm -rf /var/lib/apt/lists/*;

WORKDIR /opt/build
ADD . .

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      mkdir -p build_$TARGETARCH/src; \
      cd /opt/build/build_$TARGETARCH; \
      CC=aarch64-linux-gnu-gcc CXX=aarch64-linux-gnu-g++ cmake -DHAVE_ARM8=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all); \
      mv src/${APP} /opt/build/${APP}_arch; \
      touch /opt/build/${APP}_sse2 /opt/build/${APP}_sse41 /opt/build/${APP}_avx2; \
    else \
      if [ "$GPU" = "1" ]; then \
        export CUDACXX=/usr/local/cuda/bin/nvcc; \
        mkdir -p build_avx2/src; \
        cd /opt/build/build_avx2; \
        LIBGOMP=/usr/lib/gcc/x86_64-linux-gnu/12/; \
        cmake -GNinja -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DFORCE_STATIC_DEPS=1 \
        -DENABLE_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES="75-real;80-real;86-real;89-real;90" \
        -DOpenMP_C_FLAGS="-fopenmp -I${LIBGOMP}" -DOpenMP_C_LIB_NAMES=gomp -DOpenMP_CXX_FLAGS="-fopenmp -I${LIBGOMP}" -DOpenMP_CXX_LIB_NAMES=gomp -DOpenMP_gomp_LIBRARY=${LIBGOMP}/libgomp.a \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
        cmake --build . -j$(nproc --all); \
        mv src/${APP} /opt/build/${APP}_avx2; \
        touch /opt/build/${APP}_arch /opt/build/${APP}_sse41 /opt/build/${APP}_sse2; \
      else \
        mkdir -p build_sse2/src && mkdir -p build_sse41/src && mkdir -p build_avx2/src; \
        cd /opt/build/build_sse2; \
        cmake -DHAVE_SSE2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
        make -j $(nproc --all); \
        mv src/${APP} /opt/build/${APP}_sse2; \
        cd /opt/build/build_sse41; \
        cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
        make -j $(nproc --all); \
        mv src/${APP} /opt/build/${APP}_sse41; \
        cd /opt/build/build_avx2; \
        cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
        make -j $(nproc --all); \
        mv src/${APP} /opt/build/${APP}_avx2; \
        touch /opt/build/${APP}_arch; \
      fi; \
    fi

FROM debian:bookworm-slim
ARG TARGETARCH
ARG APP
ARG GPU

COPY --from=builder /opt/build/${APP}_arch /opt/build/${APP}_sse2 /opt/build/${APP}_sse41 /opt/build/${APP}_avx2 /usr/local/bin/
ADD util/${APP}_wrapper.sh /usr/local/bin/entrypoint

RUN apt-get update && apt-get install -y \
      gawk bash grep libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar aria2 \
    && rm -rf /var/lib/apt/lists/*; \
    if [ "$TARGETARCH" = "arm64" ]; then \
      rm -f /usr/local/bin/entrypoint; ln -s /usr/local/bin/${APP}_arch /usr/local/bin/entrypoint; \
    elif [ "$GPU" = "1" ]; then \
      rm -f /usr/local/bin/entrypoint; ln -s /usr/local/bin/${APP}_avx2 /usr/local/bin/entrypoint; \
    fi

ENTRYPOINT ["/usr/local/bin/entrypoint"]
