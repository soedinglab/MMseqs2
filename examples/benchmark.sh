#!/usr/bin/bash
set -euo pipefail

if (( $# != 5 )); then
  echo "Usage: $(basename "$0") <i:queryDb> <i:targetDb> <i:alignmentDB> <o:resultDir> <threads>" >&2
  exit 2
fi
queryDb="$1"
targetDb="$2"
alignmentDB="$3"
resultDir="$4"
threads="$5"

mkdir -p ${resultDir}/parquet/
# create TSV file
start_ns=$(date +%s%N)
../build/src/mmseqs convertalis ${queryDb} ${targetDb} ${alignmentDB} ${resultDir}/test.tsv --threads $threads --format-mode 4  --format-output query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen
status=$?
end_ns=$(date +%s%N)
elapsed_ms=$(( (end_ns - start_ns) / 1000000 ))
echo "Write to TSV time = $elapsed_ms ms" >&2

#create parquet file
start_ns=$(date +%s%N)
../build/src/mmseqs convertalis ${queryDb} ${targetDb} ${alignmentDB} ${resultDir}/parquet/test.parquet --threads $threads --format-mode 5  --format-output query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen
status=$?
end_ns=$(date +%s%N)
elapsed_ms=$(( (end_ns - start_ns) / 1000000 ))
echo "Write to parquet time = $elapsed_ms ms" >&2

# Simple DuckDB load
duckdb "benchmark.duck" <<EOF

SET threads = $threads;
.timer on

.print "Load Parquet"
SELECT count(evalue)
FROM '${resultDir}/parquet/test_*.parquet'
WHERE
    evalue BETWEEN 0.0 AND 1.0
    AND
    pident BETWEEN 25.4 AND 95.00;


.print "Load TSV"
SELECT count(evalue) 
FROM '${resultDir}/test.tsv'
WHERE
    evalue BETWEEN 0.0 AND 1.0
    AND
    pident BETWEEN 25.4 AND 95.00;
EOF

# Check if files differ
python compare_columns.py ${resultDir}/parquet ${resultDir}/test.tsv


