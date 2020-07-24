ARG NAMESPACE=
FROM debian:stable-slim as qemu-downloader
ARG NAMESPACE
RUN if [ X"$NAMESPACE" != X"" ]; then \
        apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*; \
    fi; \
    if [ X"$NAMESPACE" = X"ppc64le/" ]; then \
        wget -nv -O /usr/bin/qemu-ppc64le-static https://github.com/multiarch/qemu-user-static/releases/download/v4.2.0-4/qemu-ppc64le-static; \
        chmod +x /usr/bin/qemu-ppc64le-static; \
    fi; \
    if [ X"$NAMESPACE" = X"aarch64/" ]; then \
        wget -nv -O /usr/bin/qemu-aarch64-static https://github.com/multiarch/qemu-user-static/releases/download/v4.2.0-4/qemu-aarch64-static; \
        chmod +x /usr/bin/qemu-aarch64-static; \
    fi; \
    touch /usr/bin/dummy_copy

FROM ${NAMESPACE}debian:stable-slim as builder
ARG NAMESPACE
COPY --from=qemu-downloader /usr/bin/dummy_copy /usr/bin/qemu-aarch64-static* /usr/bin/qemu-ppc64le-static* /usr/bin/

RUN apt-get update && apt-get install -y \
    build-essential cmake xxd git zlib1g-dev libbz2-dev libatomic1 \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/mmseqs
ADD . .

RUN mkdir -p build_sse2/src && mkdir -p build_sse41/src && mkdir -p build_avx/src && mkdir -p build/src; \
    if [ X"$NAMESPACE" = X"" ]; then \
       cd /opt/mmseqs/build_sse2; \
       cmake -DHAVE_SSE2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
       make -j $(nproc --all); \
       mv src/mmseqs /opt/mmseqs/mmseqs_sse2; \
       cd /opt/mmseqs/build_sse41; \
       cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
       make -j $(nproc --all); \
       mv src/mmseqs /opt/mmseqs/mmseqs_sse41; \
       cd /opt/mmseqs/build_avx; \
       cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
       make -j $(nproc --all); \
       mv src/mmseqs /opt/mmseqs/mmseqs_avx2; \
       touch /opt/mmseqs/mmseqs_arch; \
     else \
       cd /opt/mmseqs/build; \
       cmake -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
       make -j $(nproc --all); \
       mv src/mmseqs /opt/mmseqs/mmseqs_arch; \
       touch /opt/mmseqs/mmseqs_sse2 /opt/mmseqs/mmseqs_sse42 /opt/mmseqs/mmseqs_avx2; \
     fi

FROM ${NAMESPACE}debian:stable-slim
ARG NAMESPACE
MAINTAINER Milot Mirdita <milot@mirdita.de>
COPY --from=qemu-downloader /usr/bin/dummy_copy /usr/bin/qemu-aarch64-static* /usr/bin/qemu-ppc64le-static* /usr/bin/

RUN apt-get update && apt-get install -y \
      gawk bash grep libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/mmseqs/mmseqs_arch /opt/mmseqs/mmseqs_sse2 /opt/mmseqs/mmseqs_sse41 /opt/mmseqs/mmseqs_avx2 /usr/local/bin/
ADD util/mmseqs_wrapper.sh /usr/local/bin/mmseqs
RUN if [ X"$NAMESPACE" != X"" ]; then mv -f /usr/local/bin/mmseqs_arch /usr/local/bin/mmseqs; fi

CMD ["/usr/local/bin/mmseqs"]

