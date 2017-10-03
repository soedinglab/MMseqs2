FROM alpine:latest as mmseqs-builder

RUN apk add --no-cache gcc g++ cmake musl-dev vim git ninja zlib-dev bzip2-dev

WORKDIR /opt/mmseqs
ADD . .

WORKDIR build_sse
RUN cmake -G Ninja -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
RUN ninja && ninja install

WORKDIR ../build_avx
RUN cmake -G Ninja -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
RUN ninja && ninja install

FROM alpine:latest
MAINTAINER Milot Mirdita <milot@mirdita.de>
RUN apk add --no-cache gawk bash grep libstdc++ libgomp zlib libbz2

COPY --from=mmseqs-builder /opt/mmseqs/build_sse/bin/mmseqs /usr/local/bin/mmseqs_sse42
COPY --from=mmseqs-builder /opt/mmseqs/build_avx/bin/mmseqs /usr/local/bin/mmseqs_avx2
ADD util/mmseqs_wrapper.sh /usr/local/bin/mmseqs

CMD ["mmseqs"]
