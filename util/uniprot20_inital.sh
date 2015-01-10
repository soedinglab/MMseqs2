#!/bin/sh -e

INDB="/data/hn02/share/mirdita/uniprot_full"
CLUSTERDB="/data/hn02/share/mirdita/uniprot_full_clu"
OUTDB="/data/hn02/share/mirdita/uniprot20_2014_03/uniprot20_2014_03"

TMPDIR="/tmp/hhsuite_database"
mkdir -p ${TMPDIR}

FASTADB="${TMPDIR}/hhsuite_database.$RANDOM"
#FASTADB="/data/hn02/share/mirdita/hhsuite_database.27173"
mmseqs clusteringtofastadb ${CLUSTERDB} ${INDB}_h ${INDB} ${FASTADB}

#MPIARGS="--hostfile $HOME/soeding/hostfile_cip -np 30"
MPIARGS="-np 32"
>${OUTDB}_a3m.log  2>${OUTDB}_a3m.err mpirun $MPIARGS \
 ffindex_apply_mpi ${FASTADB} ${FASTADB}.index \
 -i ${OUTDB}_a3m.ffindex -d ${OUTDB}_a3m.ffdata -r \
 -- $MMDIR/util/fasta_to_msa_a3m.sh 
>${OUTDB}_hhm.log  2>${OUTDB}_hhm.err mpirun $MPIARGS \
  ffindex_apply_mpi ${OUTDB}_a3m.ffdata ${OUTDB}_a3m.ffindex \
  -i ${OUTDB}_hhm.ffindex -d ${OUTDB}_hhm.ffdata -r \
  -- $MMDIR/util/hhmake_wrapper.sh  # & echo $! >${OUTDB}_hhm.pid

>${OUTDB}_cs219.log 2>${OUTDB}_cs219.err mpirun $MPIARGS \
  ffindex_apply_mpi ${OUTDB}_a3m.ffdata ${OUTDB}_a3m.ffindex \
  -i ${OUTDB}_cs219.ffindex -d ${OUTDB}_cs219.ffdata -r \
  -- $MMDIR/util/cstranslate_wrapper.sh # & echo $! >${OUTDB}_cs219.pid

#wait $(<${OUTDB}_hhm.pid) $(<${OUTDB}_cs219.pid)
#rm -f ${OUTDB}_hhm.pid ${OUTDB}_cs219.pid


