#!/bin/bash -ex

#clusters2ffindex uniprot_full_h uniprot_full uniprot_clustered_msa

INDIR="/home/m/mirdita/soeding"
OUTDIR="/big/m/mirdita/small"
RELEASE="uniprot20_2014_test"
#MPIARGS="--hostfile $HOME/soeding/hostfile_cip -np 30"
MPIARGS=""
mpirun $MPIARGS \
 ffindex_apply_mpi ${INDIR}/uniprot-msa ${INDIR}/uniprot-msa-head.index \
 -i ${OUTDIR}/${RELEASE}_a3m.ffindex -d ${OUTDIR}/${RELEASE}_a3m.ffdata \
 -- $MMDIR/util/fasta_to_msa_a3m.sh
mpirun $MPIARGS \
  ffindex_apply_mpi ${OUTDIR}/${RELEASE}_a3m.ffdata ${OUTDIR}/${RELEASE}_a3m.ffindex \
  -i ${OUTDIR}/${RELEASE}_hhm.ffindex -d ${OUTDIR}/${RELEASE}_hhm.ffdata \
  -- $MMDIR/util/hhmake_wrapper.sh
mpirun $MPIARGS \
  ffindex_apply_mpi ${OUTDIR}/${RELEASE}_a3m.ffdata ${OUTDIR}/${RELEASE}_a3m.ffindex \
  -i ${OUTDIR}/${RELEASE}_cs219.ffindex -d ${OUTDIR}/${RELEASE}_cs219.ffdata \
  -- $MMDIR/util/cstranslate_wrapper.sh

