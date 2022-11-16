ARG APP=mmseqs
FROM --platform=$BUILDPLATFORM debian:stable-slim as builder
ARG TARGETARCH
ARG APP

RUN dpkg --add-architecture $TARGETARCH \
    && apt-get update \
    && apt-get install -y \
      build-essential cmake xxd git \
      zlib1g-dev libbz2-dev libatomic1 \
      crossbuild-essential-$TARGETARCH zlib1g-dev:$TARGETARCH libbz2-dev:$TARGETARCH \
    && rm -rf /var/lib/apt/lists/*

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
    fi

FROM debian:stable-slim
ARG TARGETARCH
ARG APP

RUN apt-get update && apt-get install -y \
      gawk bash grep libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/build/${APP}_arch /opt/build/${APP}_sse2 /opt/build/${APP}_sse41 /opt/build/${APP}_avx2 /usr/local/bin/
ADD util/${APP}_wrapper.sh /usr/local/bin/entrypoint
RUN if [ "$TARGETARCH" = "arm64" ]; then rm -f /usr/local/bin/entrypoint; ln -s /usr/local/bin/${APP}_arch /usr/local/bin/entrypoint; fi

ENTRYPOINT ["/usr/local/bin/entrypoint"]

