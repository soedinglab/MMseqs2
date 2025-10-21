ARG TARGETARCH
ARG APP=mmseqs
ARG GPU=0

FROM --platform=$BUILDPLATFORM ghcr.io/steineggerlab/build-containers:main-sbsa AS builder-arm64

FROM --platform=$BUILDPLATFORM ghcr.io/steineggerlab/build-containers:main-x86_64 AS builder-amd64

FROM builder-${TARGETARCH} AS builder
ARG TARGETARCH
ARG APP
ARG GPU

WORKDIR /opt/build
ADD . .
RUN set -x; \
    mkdir -p build/src; \
    cd /opt/build/build; \
    if [ "$TARGETARCH" = "arm64" ]; then \
      /usr/local/bin/cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=1 -DENABLE_WERROR=1 -DHAVE_ARM8=1 \
          -DCMAKE_TOOLCHAIN_FILE=/opt/toolchain.cmake \
          -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AARCH64} -pthread" -DOpenMP_C_LIB_NAMES=omp \
          -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AARCH64} -pthread" -DOpenMP_CXX_LIB_NAMES=omp \
          -DOpenMP_C_LIB_NAMES="libomp;dl" \
          -DOpenMP_CXX_LIB_NAMES="libomp;dl" \
          -DOpenMP_libomp_LIBRARY=${LIBOMP_AARCH64}/libomp.a \
          -DOpenMP_dl_LIBRARY=dl \
          -DRust_TOOLCHAIN=stable-x86_64-unknown-linux-gnu \
          -DRust_CARGO_TARGET=aarch64-unknown-linux-gnu \
          -DCMAKE_POLICY_DEFAULT_CMP0074=NEW -DCMAKE_POLICY_DEFAULT_CMP0144=NEW \
          -DFORCE_STATIC_DEPS=1 -DENABLE_CUDA=${GPU} -DCMAKE_CUDA_ARCHITECTURES="75-real;80-real;86-real;89-real;90" ..; \
      /usr/local/bin/cmake --build . -j$(nproc --all) -v; \
    else \
      if [ -e "${LIBGCC}/libgomp.so" ]; then \
          mv -f -- "${LIBGCC}/libgomp.so" "${LIBGCC}/libgomp.so.disabled"; \
      fi; \
      /usr/local/bin/cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=1 -DENABLE_WERROR=1 -DHAVE_AVX2=1 \
          -DOpenMP_C_FLAGS="-fopenmp -I${LIBGCC} -L${LIBGCC}" -DOpenMP_C_LIB_NAMES=gomp \
          -DOpenMP_CXX_FLAGS="-fopenmp -I${LIBGCC} -L${LIBGCC}" -DOpenMP_CXX_LIB_NAMES=gomp \
          -DOpenMP_gomp_LIBRARY="${LIBGCC}/libgomp.a" \
          -DATOMIC_LIB_OVERRIDE="${LIBGCC}/libatomic.a" \
          -DCMAKE_POLICY_DEFAULT_CMP0074=NEW -DCMAKE_POLICY_DEFAULT_CMP0144=NEW \
          -DZLIB_ROOT=/deps -DBZIP2_ROOT=/deps \
          -DFORCE_STATIC_DEPS=1 -DENABLE_CUDA=${GPU} -DCMAKE_CUDA_ARCHITECTURES="75-real;80-real;86-real;89-real;90" ..; \
      /usr/local/bin/cmake --build . -j$(nproc --all) -v; \
    fi; \
    mv src/${APP} /opt/build/${APP};

FROM debian:trixie-slim
ARG TARGETARCH
ARG APP
ARG GPU

COPY --from=builder /opt/build/${APP} /usr/local/bin/

RUN apt-get update && apt-get install -y \
  gawk bash grep wget tar aria2 \
  && rm -rf /var/lib/apt/lists/*;

RUN echo "#!/bin/sh\nexec /usr/local/bin/$APP \"\${@}\"" > /usr/local/bin/entrypoint && chmod +x /usr/local/bin/entrypoint
RUN cat /usr/local/bin/entrypoint
ENTRYPOINT ["/usr/local/bin/entrypoint"]
