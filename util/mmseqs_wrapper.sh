#!/bin/bash
if $(grep -q -P '^flags.+avx2' /proc/cpuinfo); then
    exec /usr/local/bin/mmseqs_avx2 "$@"
else
    exec /usr/local/bin/mmseqs_sse42 "$@"
fi
