ARG NAMESPACE=
FROM alpine:3.7 as qemu-downloader
ARG NAMESPACE
RUN apk add --no-cache wget && \
  if [[ "$NAMESPACE" == "arm64v8/" ]]; then \
    wget -nv -O "/usr/bin/qemu-aarch64-static" https://github.com/multiarch/qemu-user-static/releases/download/v3.1.0-2/qemu-aarch64-static; \ 
  else \
    echo -e '#!/bin/sh\n"$@"\n' > "/usr/bin/qemu-aarch64-static"; \
  fi; \
  chmod +x /usr/bin/qemu-aarch64-static;

FROM ${NAMESPACE}alpine:3.7 as mmseqs-builder
ARG NAMESPACE
COPY --from=qemu-downloader /usr/bin/qemu-aarch64-static /usr/bin/qemu-aarch64-static
RUN apk add --no-cache gcc g++ cmake musl-dev vim git ninja zlib-dev bzip2-dev

WORKDIR /opt/mmseqs
ADD . .
RUN mkdir -p build_sse/bin && mkdir -p build_avx/bin && mkdir -p build_neon/bin

WORKDIR /opt/mmseqs/build_sse
RUN if [[ "$NAMESPACE" == "" ]]; then \
    cmake -G Ninja -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    ninja && ninja install; \
  fi

WORKDIR /opt/mmseqs/build_avx
RUN if [[ "$NAMESPACE" == "" ]]; then \
    cmake -G Ninja -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    ninja && ninja install; \
  fi

WORKDIR /opt/mmseqs/build_neon
RUN if [[ "$NAMESPACE" == "arm64v8/" ]]; then \
    cmake -G Ninja -DHAVE_NEON=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    ninja && ninja install; \
    touch /opt/mmseqs/build_sse/bin/mmseqs; \
    touch /opt/mmseqs/build_avx/bin/mmseqs; \
  elif [[ "$NAMESPACE" == "" ]]; then \
    touch /opt/mmseqs/build_neon/bin/mmseqs; \
  fi

FROM ${NAMESPACE}alpine:3.7
ARG NAMESPACE
MAINTAINER Milot Mirdita <milot@mirdita.de>
COPY --from=qemu-downloader /usr/bin/qemu-aarch64-static /usr/bin/qemu-aarch64-static

RUN apk add --no-cache gawk bash grep libstdc++ libgomp zlib libbz2

COPY --from=mmseqs-builder /opt/mmseqs/build_sse/bin/mmseqs /usr/local/bin/mmseqs_sse42
COPY --from=mmseqs-builder /opt/mmseqs/build_avx/bin/mmseqs /usr/local/bin/mmseqs_avx2
COPY --from=mmseqs-builder /opt/mmseqs/build_neon/bin/mmseqs /usr/local/bin/mmseqs_neon
ADD util/mmseqs_wrapper.sh /usr/local/bin/mmseqs

RUN if [[ "$NAMESPACE" == "arm64v8/" ]]; then mv -f /usr/local/bin/mmseqs_neon /usr/local/bin/mmseqs; fi

CMD ["/usr/local/bin/mmseqs"]

