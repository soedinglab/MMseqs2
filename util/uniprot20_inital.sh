#!/bin/sh -e

INDB="/home/m/mirdita/soeding/mmseqs/scopff"
CLUSTERDB="/home/m/mirdita/soeding/mmseqs/test"
OUTDB="/big/m/mirdita/small/scopff_test"

TMPDIR="/tmp/hhsuite_database"
mkdir -p ${TMPDIR}

FASTADB="${TMPDIR}/hhsuite_database.$RANDOM"
mmseqs clusteringtofastadb ${CLUSTERDB} ${INDB}_h ${INDB} ${FASTADB}

#MPIARGS="--hostfile $HOME/soeding/hostfile_cip -np 30"
MPIARGS=""
mpirun $MPIARGS \
 ffindex_apply_mpi ${FASTADB} ${FASTADB}.index \
 -i ${OUTDB}_a3m.ffindex -d ${OUTDB}_a3m.ffdata \
 -- $MMDIR/util/fasta_to_msa_a3m.sh
mpirun $MPIARGS \
  ffindex_apply_mpi ${OUTDB}_a3m.ffdata ${OUTDB}_a3m.ffindex \
  -i ${OUTDB}_hhm.ffindex -d ${OUTDB}_hhm.ffdata \
  -- $MMDIR/util/hhmake_wrapper.sh
mpirun $MPIARGS \
  ffindex_apply_mpi ${OUTDB}_a3m.ffdata ${OUTDB}_a3m.ffindex \
  -i ${OUTDB}_cs219.ffindex -d ${OUTDB}_cs219.ffdata \
  -- $MMDIR/util/cstranslate_wrapper.sh

