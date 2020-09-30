#!/bin/sh
FLAGS="$(grep -m 1 '^flags' /proc/cpuinfo)"
case "${FLAGS}" in
  *avx2*)
    exec /usr/local/bin/mmseqs_avx2 "$@"
    ;;
  *sse4_1*)
    exec /usr/local/bin/mmseqs_sse41 "$@"
    ;;
  *)
    exec /usr/local/bin/mmseqs_sse2 "$@"
    ;;
esac
