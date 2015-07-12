#!/bin/sh -ex
#BSUB -q mpi
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -W 48:00
#BSUB -m hh
#BSUB -a openmpi 
#BSUB -n 192 
#BSUB -R haswell 
#BSUB -R cbscratch 
#BSUB -R "span[ptile=16]"
INDB="/cbscratch/martin/unref50/uniref50"
CLUSTERDB="/cbscratch/martin/unref50/uniref50_CLU"
OUTDB="/cbscratch/martin/unref50/uniref_hhblits"

TMPDIR="/cbscratch/martin/hhsuite_database"
mkdir -p ${TMPDIR}
MPIRUNNER=mpirun
FASTADB="${TMPDIR}/hhsuite_database.$RANDOM"
#FASTADB="/cbscratch/martin/unref50/uniref50"
mmseqs clustertofastadb ${CLUSTERDB} ${INDB}_h ${INDB} ${FASTADB}

#MPIARGS="--hostfile $HOME/soeding/hostfile_cip -np 30"
#MPIARGS="-np 32"
MPIARGS=""
>${OUTDB}_a3m.log  2>${OUTDB}_a3m.err $MPIRUNNER $MPIARGS \
 ffindex_apply_mpi ${FASTADB} ${FASTADB}.index \
 -i ${OUTDB}_a3m.ffindex -d ${OUTDB}_a3m.ffdata \
 -- $MMDIR/util/fasta_to_msa_a3m.sh 

>${OUTDB}_hhm.log  2>${OUTDB}_hhm.err $MPIRUNNER $MPIARGS \
  ffindex_apply_mpi ${OUTDB}_a3m.ffdata ${OUTDB}_a3m.ffindex \
  -i ${OUTDB}_hhm.ffindex -d ${OUTDB}_hhm.ffdata \
  -- $MMDIR/util/hhmake_wrapper.sh  # & echo $! >${OUTDB}_hhm.pid

>${OUTDB}_cs219.log 2>${OUTDB}_cs219.err $MPIRUNNER $MPIARGS \
  ffindex_apply_mpi ${OUTDB}_a3m.ffdata ${OUTDB}_a3m.ffindex \
  -i ${OUTDB}_cs219.ffindex -d ${OUTDB}_cs219.ffdata \
  -- $MMDIR/util/cstranslate_wrapper.sh # & echo $! >${OUTDB}_cs219.pid

#wait $(<${OUTDB}_hhm.pid) $(<${OUTDB}_cs219.pid)
#rm -f ${OUTDB}_hhm.pid ${OUTDB}_cs219.pid

mkdir -p ordered
OUTNAME=$(basename ${OUTDB})
# order by lenght for speed up computation (viterbi)
sort -n -k3 -r ${OUTDB}_cs219.ffindex | grep -v "[[:space:]]1$"  | awk '{print $1}' > ${OUTDB}.length

ffindex_order ${OUTDB}.length ${OUTDB}_a3m.ffdata ${OUTDB}_a3m.ffindex ordered/${OUTNAME}_a3m.ffdata ordered/${OUTNAME}_a3m.ffindex
ffindex_order ${OUTDB}.length ${OUTDB}_cs219.ffdata ${OUTDB}_cs219.ffindex ordered/${OUTNAME}_cs219.ffdata ordered/${OUTNAME}_cs219.ffindex
ffindex_order ${OUTDB}.length ${OUTDB}_hhm.ffdata ${OUTDB}_hhm.ffindex ordered/${OUTNAME}_hhm.ffdata ordered/${OUTNAME}_hhm.ffindex

# remove hhms with size of 1 
awk '{if($3 >= 20) {print $0}}' ${OUTDB}_hhm.ffindex > valid_hhm
ffindex_order valid_hhm ordered/${OUTNAME}_hhm.ffdata ordered/${OUTNAME}_hhm.ffindex ordered/${OUTNAME}_hhm2.ffdata ordered/${OUTNAME}_hhm2.ffindex
ffindex_order ${OUTDB}.length ordered/${OUTNAME}_hhm2.ffdata ordered/${OUTNAME}_hhm2.ffindex ordered/${OUTNAME}_hhm.ffdata ordered/${OUTNAME}_hhm.ffindex

# support for old HH-suite (2.0.16)
# generate files for legacy hh-suite compatibility
mmseqs legacycs219 ${OUTDB}_cs219 ${OUTDB}_a3m ${OUTDB}

ln -s ${OUTDB}_a3m.ffdata ${OUTDB}_a3m_db
ln -s ${OUTDB}_hhm.ffdata ${OUTDB}_hhm_db

# change file endings 
awk '{print $1".a3m\t"$2"\t"$3}' ${OUTDB}_a3m.ffindex > ${OUTDB}_a3m_db.index
awk '{print $1".hhm\t"$2"\t"$3}' ${OUTDB}_hhm.ffindex > ${OUTDB}_hhm_db.index
