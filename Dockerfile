ARG NAMESPACE=
FROM debian:stable-slim as qemu-downloader
ARG NAMESPACE
RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*
RUN if [ X"$NAMESPACE" = X"arm64v8/" ]; then \
      wget -nv -O "/usr/bin/qemu-aarch64-static" https://github.com/multiarch/qemu-user-static/releases/download/v3.1.0-2/qemu-aarch64-static; \
    else \
      echo -e '#!/bin/sh\n"$@"\n' > "/usr/bin/qemu-aarch64-static"; \
    fi; \
    chmod +x /usr/bin/qemu-aarch64-static;

FROM ${NAMESPACE}debian:stable-slim as mmseqs-builder
ARG NAMESPACE
COPY --from=qemu-downloader /usr/bin/qemu-aarch64-static /usr/bin/qemu-aarch64-static

RUN apt-get update && apt-get install -y \
    build-essential cmake xxd git zlib1g-dev libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/mmseqs
ADD . .
RUN mkdir -p build_sse/bin && mkdir -p build_avx/bin && mkdir -p build_neon/bin

WORKDIR /opt/mmseqs/build_sse
RUN if [ X"$NAMESPACE" = X"" ]; then \
      cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all) && make install; \
    fi

WORKDIR /opt/mmseqs/build_avx
RUN if [ X"$NAMESPACE" = X"" ]; then \
      cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all) && make install; \
    fi

WORKDIR /opt/mmseqs/build_neon
RUN if [ X"$NAMESPACE" = X"arm64v8/" ]; then \
      cmake  -DHAVE_NEON=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all) && make install; \
      touch /opt/mmseqs/build_sse/bin/mmseqs; \
      touch /opt/mmseqs/build_avx/bin/mmseqs; \
    else \
      touch /opt/mmseqs/build_neon/bin/mmseqs; \
    fi

FROM ${NAMESPACE}debian:stable-slim
ARG NAMESPACE
MAINTAINER Milot Mirdita <milot@mirdita.de>
COPY --from=qemu-downloader /usr/bin/qemu-aarch64-static /usr/bin/qemu-aarch64-static

RUN apt-get update && apt-get install -y \
      gawk bash grep libstdc++6 libgomp1 zlib1g libbz2-1.0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=mmseqs-builder /opt/mmseqs/build_sse/bin/mmseqs /usr/local/bin/mmseqs_sse42
COPY --from=mmseqs-builder /opt/mmseqs/build_avx/bin/mmseqs /usr/local/bin/mmseqs_avx2
COPY --from=mmseqs-builder /opt/mmseqs/build_neon/bin/mmseqs /usr/local/bin/mmseqs_neon
ADD util/mmseqs_wrapper.sh /usr/local/bin/mmseqs

RUN if [ X"$NAMESPACE" = X"arm64v8/" ]; then mv -f /usr/local/bin/mmseqs_neon /usr/local/bin/mmseqs; fi

CMD ["/usr/local/bin/mmseqs"]

