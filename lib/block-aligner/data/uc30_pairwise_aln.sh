mmseqs align seqdb seqdb uc30_clu uc30_aln --comp-bias-corr 0 -a --threads 128
mmseqs convertalis seqdb seqdb uc30_aln uc30_aln.m8 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qseq,tseq --threads 16
awk '$1!=prev && $3 < 0.4 {print; prev=$1}' uc30_aln.m8 | shuf | head -n 1000 > uc30_30_40.m8
awk '$1!=prev && $3 > 0.4 && $3 < 0.5 {print; prev=$1}' uc30_aln.m8 | shuf | head -n 1000 > uc30_40_50.m8
awk '$1!=prev && $3 > 0.5 && $3 < 0.6 {print; prev=$1}' uc30_aln.m8 | shuf | head -n 1000 > uc30_50_60.m8
awk '$1!=prev && $3 > 0.6 && $3 < 0.7 {print; prev=$1}' uc30_aln.m8 | shuf | head -n 1000 > uc30_60_70.m8
awk '$1!=prev && $3 > 0.7 && $3 < 0.8 {print; prev=$1}' uc30_aln.m8 | shuf | head -n 1000 > uc30_70_80.m8
awk '$1!=prev && $3 > 0.8 && $3 < 0.9 {print; prev=$1}' uc30_aln.m8 | shuf | head -n 1000 > uc30_80_90.m8
awk '$1!=prev && $3 > 0.9 && $3 < 1.0 {print; prev=$1}' uc30_aln.m8 | shuf | head -n 1000 > uc30_90_100.m8
