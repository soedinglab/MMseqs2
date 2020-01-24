#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
    [ ! -f "$1" ]
}

hasCommand () {
    command -v "$1" >/dev/null 2>&1
}

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
    fail "No download tool found in PATH. Please install either aria2c, curl or wget."
fi

downloadFile() {
    URL="$1"
    OUTPUT="$2"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            aria2c --max-connection-per-server="$THREAD_NUM" --allow-overwrite=true -o "$OUTPUT" "$URL" && return 0
            ;;
        CURL)
            curl -o "$OUTPUT" "$URL" && return 0
            ;;
        WGET)
            wget -O "$OUTPUT" "$URL" && return 0
            ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

# check number of input variables
[ "$#" -ne 3 ] && echo "Please provide <selection> <outDB> <tmp>" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";

SELECTION="$1"
OUTDB="$2"
TMP_PATH="$3"

INPUT_TYPE=""
case "${SELECTION}" in
    "UniRef100")
        if notExists "${TMP_PATH}/db.fasta.gz"; then
          downloadFile "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.release_note" "${OUTDB}.version"
          downloadFile "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz" "${TMP_PATH}/db.fasta.gz"
          INPUT_TYPE="AA"
        fi
    ;;
    "UniRef90")
        if notExists "${TMP_PATH}/db.fasta.gz"; then
          downloadFile "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.release_note" "${OUTDB}.version"
          downloadFile "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz" "${TMP_PATH}/db.fasta.gz"
          INPUT_TYPE="AA"
        fi
    ;;
    "UniRef50")
        if notExists "${TMP_PATH}/db.fasta.gz"; then
          downloadFile "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.release_note" "${OUTDB}.version"
          downloadFile "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz" "${TMP_PATH}/db.fasta.gz"
          INPUT_TYPE="AA"
        fi
    ;;
    "UniProtKB")
        if notExists "${TMP_PATH}/db1.fasta.gz" || notExists "${TMP_PATH}/db2.fasta.gz"; then
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt" "${OUTDB}.version"
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" "${TMP_PATH}/db1.fasta.gz"
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz" "${TMP_PATH}/db2.fasta.gz"
          INPUT_TYPE="AA_2"
        fi
    ;;
    "UniProtKB/TrEMBL")
        if notExists "${TMP_PATH}/db.fasta.gz"; then
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt" "${OUTDB}.version"
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz" "${TMP_PATH}/db.fasta.gz"
          INPUT_TYPE="AA"
        fi
    ;;
    "UniProtKB/Swiss-Prot")
        if notExists "${TMP_PATH}/db.fasta.gz"; then
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt" "${OUTDB}.version"
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" "${TMP_PATH}/db.fasta.gz"
          INPUT_TYPE="AA"
        fi
    ;;
    "NR")
        if notExists "${TMP_PATH}/db.fasta.gz"; then
          date "+%s" >> "${OUTDB}.version"
          downloadFile "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz" "${TMP_PATH}/db.fasta.gz"
          INPUT_TYPE="AA"
        fi
    ;;
    "NT")
        if notExists "${TMP_PATH}/db.fasta.gz"; then
          date "+%s" >> "${OUTDB}.version"
          downloadFile "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz" "${TMP_PATH}/db.fasta.gz"
          INPUT_TYPE="AA"
        fi
    ;;
    "PDB")
        if notExists "${TMP_PATH}/db.fasta.gz"; then
          date "+%s" >> "${OUTDB}.version"
          downloadFile "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz" "${TMP_PATH}/db.fasta.gz"
          INPUT_TYPE="AA"
        fi
    ;;
    "PDB70")
        if notExists "${TMP_PATH}/msa.index"; then
          date "+%s" >> "${OUTDB}.version"
          downloadFile "http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz" "${TMP_PATH}/pdb70.tar.gz"
          tar -xOzf "${TMP_PATH}/pdb70.tar.gz" pdb70_a3m.ffdata | tr -d '\000' | awk -v outfile="${TMP_PATH}/msa" 'function writeEntry() { printf "%s\0", data >> outfile; size = length(data) + 1; data=""; print id"\t"offset"\t"size >> outindex; offset = offset + size; } BEGIN { data = ""; offset = 0; id = 1; if(length(outfile) == 0) { outfile="output"; } outindex = outfile".index"; printf("") > outfile; printf("") > outindex; } /^>ss_/ { inss = 1; entry = 0; next; } inss == 1 { inss = 0; next; } /^>/ && entry == 0 { if (id > 1) { writeEntry(); } id = id + 1; data = ">"substr($1, 2)"\n"; entry = entry + 1; next; } entry > 0 { data = data""$0"\n"; entry = entry + 1; next; } END { writeEntry(); close(outfile); close(outfile".index"); }'
          rm -f "${TMP_PATH}/pdb70.tar.gz"
          INPUT_TYPE="A3M"
        fi
    ;;
    "Pfam-A.full")
        if notExists "${TMP_PATH}/db.msa.gz"; then
          downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz" "${OUTDB}.version"
          downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz" "${TMP_PATH}/db.msa.gz"
          INPUT_TYPE="MSA"
        fi
    ;;
    "Pfam-A.seed")
        if notExists "${TMP_PATH}/db.msa.gz"; then
          downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz" "${OUTDB}.version"
          downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz" "${TMP_PATH}/db.msa.gz"
          INPUT_TYPE="MSA"
        fi
    ;;
    "Resfinder")
        if notExists "${TMP_PATH}/download.done"; then
          downloadFile "https://api.bitbucket.org/2.0/repositories/genomicepidemiology/resfinder_db/commit/master?fields=hash,date" "${OUTDB}.version"
          downloadFile "https://bitbucket.org/genomicepidemiology/resfinder_db/get/master.tar.gz" "${TMP_PATH}/master.tar.gz"
          tar -C "${TMP_PATH}" --strip-components=1 -xzvf "${TMP_PATH}/master.tar.gz" "*.fsa"
          touch "${TMP_PATH}/download.done"
          rm -f "${TMP_PATH}/master.tar.gz"
          INPUT_TYPE="FSA"
        fi
    ;;

esac

if notExists "${OUTDB}.dbtype"; then
case "${INPUT_TYPE}" in
    "AA")
        # shellcheck disable=SC2086
        "${MMSEQS}" createdb "${TMP_PATH}/db.fasta.gz" "${OUTDB}" ${COMP_PAR} \
            || fail "createdb died"
    ;;
    "AA_2")
        # shellcheck disable=SC2086
        "${MMSEQS}" createdb "${TMP_PATH}/db1.fasta.gz" "${TMP_PATH}/db2.fasta.gz" "${OUTDB}" ${COMP_PAR} \
            || fail "createdb died"
    ;;
    "FSA")
        # shellcheck disable=SC2086
        "${MMSEQS}" createdb "${TMP_PATH}/"*.fsa "${OUTDB}" ${COMP_PAR} \
            || fail "createdb died"
    ;;
    "A3M")
        # shellcheck disable=SC2086
        "${MMSEQS}" msa2profile "${TMP_PATH}/msa" "${OUTDB}" --match-mode 1 --match-ratio 0.5 --msa-type 1 ${THREADS_PAR} \
            || fail "msa2profile died"
    ;;
    "MSA")
        # shellcheck disable=SC2086
        "${MMSEQS}" convertmsa "${TMP_PATH}/db.msa.gz" "${TMP_PATH}/msa" ${VERB_PAR} \
            || fail "convertmsa died"
        # shellcheck disable=SC2086
        "${MMSEQS}" msa2profile "${TMP_PATH}/msa" "${OUTDB}" --match-mode 1 --match-ratio 0.5 ${THREADS_PAR} \
            || fail "msa2profile died"
    ;;
esac
fi
if [ -n "${TAXONOMY}" ] && notExists "${OUTDB}_mapping"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" prefixid "${OUTDB}_h" "${TMP_PATH}/header_pref.tsv" --tsv ${THREADS_PAR} \
        || fail "prefixid died"
    awk '{match($0, /(OX|TaxID)=[0-9]+ /); print $1"\t"substr($0, RSTART+3, RLENGTH-4); }' "${TMP_PATH}/header_pref.tsv" \
        | LC_ALL=C sort -n > "${OUTDB}_mapping"
    # shellcheck disable=SC2086
    "${MMSEQS}" createtaxdb "${OUTDB}" "${TMP_PATH}/taxonomy" ${THREADS_PAR} \
        || fail "createtaxdb died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    rm -f "${TMP_PATH}/db.fasta.gz"
    rm -f "${TMP_PATH}/db1.fasta.gz" "${TMP_PATH}/db2.fasta.gz"
    rm -f "${TMP_PATH}/db.msa.gz"
    rm -f "${TMP_PATH}/msa" "${TMP_PATH}/msa.index"
    rm -f "${TMP_PATH}/header_pref.tsv"
    rm -f "${TMP_PATH}/download.sh"
fi
