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

ARR=""
push_back() {
    # shellcheck disable=SC1003
    CURR="$(printf '%s' "$1" | awk '{ gsub(/'\''/, "'\''\\'\'''\''"); print; }')"
    if [ -z "$ARR" ]; then
        ARR=''\'$CURR\'''
    else
        ARR=$ARR' '\'$CURR\'''
    fi
}

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
    fail "No download tool found in PATH. Please install aria2c, curl or wget."
fi

downloadFile() {
    URL="$1"
    OUTPUT="$2"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            FILENAME=$(basename "${OUTPUT}")
            DIR=$(dirname "${OUTPUT}")
            aria2c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -o "$FILENAME" -d "$DIR" "$URL" && return 0
            ;;
        CURL)
            curl -L -o "$OUTPUT" "$URL" && return 0
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
        if notExists "${TMP_PATH}/uniref100.fasta.gz"; then
            downloadFile "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.release_note" "${TMP_PATH}/version"
            downloadFile "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz" "${TMP_PATH}/uniref100.fasta.gz"
        fi
        push_back "${TMP_PATH}/uniref100.fasta.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "UniRef90")
        if notExists "${TMP_PATH}/uniref90.fasta.gz"; then
            downloadFile "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.release_note" "${TMP_PATH}/version"
            downloadFile "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz" "${TMP_PATH}/uniref90.fasta.gz"
        fi
        push_back "${TMP_PATH}/uniref90.fasta.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "UniRef50")
        if notExists "${TMP_PATH}/uniref50.fasta.gz"; then
            downloadFile "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.release_note" "${TMP_PATH}/version"
            downloadFile "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz" "${TMP_PATH}/uniref50.fasta.gz"
        fi
        push_back "${TMP_PATH}/uniref50.fasta.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "UniProtKB")
        if notExists "${TMP_PATH}/uniprot_sprot.fasta.gz" || notExists "${TMP_PATH}/uniprot_trembl.fasta.gz"; then
            downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt" "${TMP_PATH}/version"
            downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" "${TMP_PATH}/uniprot_sprot.fasta.gz"
            downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz" "${TMP_PATH}/uniprot_trembl.fasta.gz"
        fi
        push_back "${TMP_PATH}/uniprot_sprot.fasta.gz"
        push_back "${TMP_PATH}/uniprot_trembl.fasta.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "UniProtKB/TrEMBL")
        if notExists "${TMP_PATH}/uniprot_trembl.fasta.gz"; then
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt" "${TMP_PATH}/version"
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz" "${TMP_PATH}/uniprot_trembl.fasta.gz"
        fi
        push_back "${TMP_PATH}/uniprot_trembl.fasta.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "UniProtKB/Swiss-Prot")
        if notExists "${TMP_PATH}/uniprot_sprot.fasta.gz"; then
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/reldate.txt" "${TMP_PATH}/version"
          downloadFile "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" "${TMP_PATH}/uniprot_sprot.fasta.gz"
        fi
        push_back "${TMP_PATH}/uniprot_sprot.fasta.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "NR")
        if notExists "${TMP_PATH}/nr.gz"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz" "${TMP_PATH}/nr.gz"
            downloadFile "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz" "${TMP_PATH}/prot.accession2taxid.gz"
            gunzip "${TMP_PATH}/prot.accession2taxid.gz"
            downloadFile "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz" "${TMP_PATH}/pdb.accession2taxid.gz"
            gunzip "${TMP_PATH}/pdb.accession2taxid.gz"
        fi
        push_back "${TMP_PATH}/nr.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "NT")
        if notExists "${TMP_PATH}/nt.gz"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz" "${TMP_PATH}/nt.gz"
        fi
        push_back "${TMP_PATH}/nt.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "GTDB")
        if notExists "${TMP_PATH}/download.done"; then
            downloadFile "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/VERSION.txt" "${TMP_PATH}/version"
            downloadFile "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz" "${TMP_PATH}/gtdb.tar.gz"
            downloadFile "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_taxonomy.tsv" "${TMP_PATH}/bac120_taxonomy.tsv"
            downloadFile "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_taxonomy.tsv" "${TMP_PATH}/ar53_taxonomy.tsv"
            touch "${TMP_PATH}/download.done"
        fi
        INPUT_TYPE="GTDB"
    ;;
    "PDB")
        if notExists "${TMP_PATH}/pdb_seqres.txt.gz"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz" "${TMP_PATH}/pdb_seqres.txt.gz"
        fi
        push_back "${TMP_PATH}/pdb_seqres.txt.gz"
        INPUT_TYPE="FASTA_LIST"
    ;;
    "PDB70")
        if notExists "${TMP_PATH}/msa.index"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz" "${TMP_PATH}/pdb70.tar.gz"
            tar -xOzf "${TMP_PATH}/pdb70.tar.gz" pdb70_a3m.ffdata | tr -d '\000' | awk -v outfile="${TMP_PATH}/msa" 'function writeEntry() { printf "%s\0", data >> outfile; size = length(data) + 1; data=""; print id"\t"offset"\t"size >> outindex; offset = offset + size; } BEGIN { data = ""; offset = 0; id = 1; if(length(outfile) == 0) { outfile="output"; } outindex = outfile".index"; printf("") > outfile; printf("") > outindex; printf("%c%c%c%c",11,0,0,0) > outfile".dbtype"; } /^>ss_/ { inss = 1; entry = 0; next; } inss == 1 { inss = 0; next; } /^>/ && entry == 0 { if (id > 1) { writeEntry(); } id = id + 1; data = ">"substr($1, 2)"\n"; entry = entry + 1; next; } entry > 0 { data = data""$0"\n"; entry = entry + 1; next; } END { writeEntry(); close(outfile); close(outfile".index"); }'
            rm -f "${TMP_PATH}/pdb70.tar.gz"
        fi
        INPUT_TYPE="A3M"
    ;;
    "Pfam-A.full")
        if notExists "${TMP_PATH}/db.msa.gz"; then
            downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz" "${TMP_PATH}/version"
            downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz" "${TMP_PATH}/db.msa.gz"
        fi
        INPUT_TYPE="STOCKHOLM_MSA"
    ;;
    "Pfam-A.seed")
        if notExists "${TMP_PATH}/db.msa.gz"; then
            downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz" "${TMP_PATH}/version"
            downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz" "${TMP_PATH}/db.msa.gz"
        fi
        INPUT_TYPE="STOCKHOLM_MSA"
    ;;
    "Pfam-B")
        if notExists "${TMP_PATH}/msa.tar.gz"; then
            downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz" "${TMP_PATH}/version"
            downloadFile "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-B.tgz" "${TMP_PATH}/msa.tar.gz"
        fi
        INPUT_TYPE="FASTA_MSA"
    ;;
    "eggNOG")
        if notExists "${TMP_PATH}/download.done"; then
            date "+%s" > "${TMP_PATH}/version"
            downloadFile "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2/2_raw_algs.tar" "${TMP_PATH}/bacteria.tar"
            downloadFile "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2157/2157_raw_algs.tar" "${TMP_PATH}/archea.tar"
            downloadFile "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2759/2759_raw_algs.tar" "${TMP_PATH}/eukaryota.tar"
            downloadFile "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/10239/10239_raw_algs.tar" "${TMP_PATH}/viruses.tar"
            touch "${TMP_PATH}/download.done"
        fi
        INPUT_TYPE="eggNOG"
        push_back "${TMP_PATH}/bacteria.tar"
        push_back "${TMP_PATH}/archea.tar"
        push_back "${TMP_PATH}/eukaryota.tar"
        push_back "${TMP_PATH}/viruses.tar"
        TAR2DB_INCLUDE='\.raw_alg\.faa\.gz$'
        SED_FIX_LOOKUP='s|\.raw_alg\.faa\.gz||g'
    ;;
    "VOGDB")
        if notExists "${TMP_PATH}/download.done"; then
            downloadFile "http://fileshare.csb.univie.ac.at/vog/latest/release.txt" "${TMP_PATH}/version"
            downloadFile "http://fileshare.csb.univie.ac.at/vog/latest/vog.raw_algs.tar.gz" "${TMP_PATH}/vog.tar.gz"
            touch "${TMP_PATH}/download.done"
        fi
        INPUT_TYPE="eggNOG"
        push_back "${TMP_PATH}/vog.tar.gz"
        TAR2DB_INCLUDE='\.msa$'
        SED_FIX_LOOKUP='s|\.msa||g'
    ;;
    "CDD")
        if notExists "${TMP_PATH}/msa.msa.gz"; then
            downloadFile "https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.info" "${TMP_PATH}/version"
            downloadFile "https://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz" "${TMP_PATH}/msa.tar.gz"
        fi
        INPUT_TYPE="FASTA_MSA"
        SED_FIX_LOOKUP='s|\.FASTA||g'
        FASTA_MSA_MSA2PROFILE_PAR="--skip-query"
    ;;
    "Resfinder")
        if notExists "${TMP_PATH}/download.done"; then
            downloadFile "https://api.bitbucket.org/2.0/repositories/genomicepidemiology/resfinder_db/commit/master?fields=hash,date" "${TMP_PATH}/version"
            downloadFile "https://bitbucket.org/genomicepidemiology/resfinder_db/get/master.tar.gz" "${TMP_PATH}/master.tar.gz"
            # avoid tar wildcard extraction as it's not available in busybox tar (windows, biocontainer)
            mkdir -p "${TMP_PATH}/fsa"
            tar -C "${TMP_PATH}/fsa" --strip-components=1 -xzvf "${TMP_PATH}/master.tar.gz"
            mv -f -- "${TMP_PATH}/fsa/"*.fsa "${TMP_PATH}"
            rm -rf -- "${TMP_PATH}/master.tar.gz" "${TMP_PATH}/fsa"
            touch "${TMP_PATH}/download.done"
        fi
        INPUT_TYPE="FSA"
    ;;
    "dbCAN2")
        if notExists "${TMP_PATH}/download.done"; then
            downloadFile "http://bcb.unl.edu/dbCAN2/download/dbCAN-fam-aln-V9.tar.gz" "${TMP_PATH}/msa.tar.gz"
            printf "9 %s\n" "$(date "+%s")" > "${TMP_PATH}/version"
            touch "${TMP_PATH}/download.done"
        fi
        INPUT_TYPE="FASTA_MSA"
        SED_FIX_LOOKUP='s|\.aln||g'
    ;;
    "SILVA")
       if notExists "${TMP_PATH}/download.done"; then
         downloadFile "https://ftp.arb-silva.de/release_${SILVA_REL}/Exports/README.txt" "${TMP_PATH}/version"
         downloadFile "https://ftp.arb-silva.de/release_${SILVA_REL}/Exports/SILVA_${SILVA_REL}_SSURef_NR99_tax_silva.fasta.gz" "${TMP_PATH}/silva.fasta.gz"
         downloadFile "https://ftp.arb-silva.de/release_${SILVA_REL}/Exports/taxonomy/tax_slv_ssu_${SILVA_REL}.txt.gz" "${TMP_PATH}/silva_tax.txt.gz"
         gunzip "${TMP_PATH}/silva_tax.txt.gz"
         downloadFile "https://ftp.arb-silva.de/release_${SILVA_REL}/Exports/taxonomy/tax_slv_ssu_${SILVA_REL}.acc_taxid.gz" "${TMP_PATH}/silva.acc_taxid.gz"
         gunzip "${TMP_PATH}/silva.acc_taxid.gz"
         touch "${TMP_PATH}/download.done"
      fi
      push_back "${TMP_PATH}/silva.fasta.gz"
      INPUT_TYPE="FASTA_LIST"
    ;;
    "Kalamari")
        if notExists "${TMP_PATH}/kalamari.tsv"; then
            printf "3.7 %s\n" "$(date "+%s")" > "${TMP_PATH}/version"
            downloadFile "https://raw.githubusercontent.com/lskatz/Kalamari/18d71da740546ba4a5117682e1ae2a037379afe0/src/Kalamari_v3.7.tsv" "${TMP_PATH}/kalamari.tsv"
        fi
        ACCESSIONS=""
        # shellcheck disable=SC2034
        while IFS="$(printf '\t')" read -r NAME ACCESSION TAXID; do
            if [ "$NAME" = "scientificName" ]; then
                continue
            fi
            case "${ACCESSION}" in XXX*)
                continue
            esac
            if [ -z "$ACCESSIONS" ]; then
                ACCESSIONS="$ACCESSION"
            else
                ACCESSIONS="$ACCESSIONS,$ACCESSION"
            fi
        done < "${TMP_PATH}/kalamari.tsv"
        if notExists "${TMP_PATH}/kalamari.fasta"; then
            # Reset download strategy to not use aria2c for NCBI
            STRATEGY=""
            if hasCommand curl; then STRATEGY="$STRATEGY CURL"; fi
            if hasCommand wget; then STRATEGY="$STRATEGY WGET"; fi
            if [ "$STRATEGY" = "" ]; then
                fail "No download tool found in PATH. Please install aria2c, curl or wget."
            fi
            downloadFile "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACCESSIONS}&rettype=fasta&retmode=txt" "${TMP_PATH}/kalamari.fasta.tmp"
            awk '/<!DOCTYPE/ { exit 1; } ' "${TMP_PATH}/kalamari.fasta.tmp" || fail "Could not download genomes from NCBI. Please try again later."
            awk -F '[\t>.]' 'NR == FNR { f[$2]=$NF; next; } /^>/{ print $0" TaxID="f[$2]" "; next; } { print; }' "${TMP_PATH}/kalamari.tsv" "${TMP_PATH}/kalamari.fasta.tmp" > "${TMP_PATH}/kalamari.fasta"
            rm -f "${TMP_PATH}/kalamari.fasta.tmp"
        fi
        push_back "${TMP_PATH}/kalamari.fasta"
        INPUT_TYPE="FASTA_LIST"
    ;;
esac

if notExists "${OUTDB}.dbtype"; then
case "${INPUT_TYPE}" in
    "FASTA_LIST")
        eval "set -- $ARR"
        # shellcheck disable=SC2086
        "${MMSEQS}" createdb "${@}" "${OUTDB}" ${COMP_PAR} \
            || fail "createdb died"
        for i in "${@}"; do
            rm -f -- "$i"
        done
    ;;
    "FSA")
        # shellcheck disable=SC2086
        "${MMSEQS}" createdb "${TMP_PATH}/"*.fsa "${OUTDB}" ${COMP_PAR} \
            || fail "createdb died"
        rm -f -- "${TMP_PATH}/"*.fsa
    ;;
    "A3M")
        # shellcheck disable=SC2086
        "${MMSEQS}" msa2profile "${TMP_PATH}/msa" "${OUTDB}" --match-mode 1 --match-ratio 0.5 --msa-type 1 ${THREADS_PAR} \
            || fail "msa2profile died"
    ;;
    "STOCKHOLM_MSA")
        # shellcheck disable=SC2086
        "${MMSEQS}" convertmsa "${TMP_PATH}/db.msa.gz" "${TMP_PATH}/msa" ${VERB_PAR} \
            || fail "convertmsa died"
        rm -f "${TMP_PATH}/db.msa.gz"
        # shellcheck disable=SC2086
        "${MMSEQS}" msa2profile "${TMP_PATH}/msa" "${OUTDB}" --match-mode 1 --match-ratio 0.5 ${THREADS_PAR} \
            || fail "msa2profile died"
        if [ -n "${REMOVE_TMP}" ]; then
            # shellcheck disable=SC2086
            "${MMSEQS}" rmdb "${TMP_PATH}/msa" ${VERB_PAR} \
                || fail "rmdb died"
        fi
    ;;
    "FASTA_MSA")
        # shellcheck disable=SC2086
        "${MMSEQS}" tar2db "${TMP_PATH}/msa.tar.gz" "${TMP_PATH}/msa" --output-dbtype 11 ${THREADS_PAR}  \
            || fail "tar2db died"
        if [ -n "${SED_FIX_LOOKUP}" ]; then
            sed "${SED_FIX_LOOKUP}" "${TMP_PATH}/msa.lookup" > "${TMP_PATH}/msa.lookup_tmp"
            mv -f -- "${TMP_PATH}/msa.lookup_tmp" "${TMP_PATH}/msa.lookup"
        fi
        rm -f -- "${TMP_PATH}/msa.tar.gz"
        # shellcheck disable=SC2086
        "${MMSEQS}" msa2profile "${TMP_PATH}/msa" "${OUTDB}" --match-mode 1 --match-ratio 0.5 ${FASTA_MSA_MSA2PROFILE_PAR} ${THREADS_PAR} \
            || fail "msa2profile died"
        if [ -n "${REMOVE_TMP}" ]; then
            # shellcheck disable=SC2086
            "${MMSEQS}" rmdb "${TMP_PATH}/msa" ${VERB_PAR} \
                || fail "rmdb died"
        fi
    ;;
    "eggNOG")
        eval "set -- $ARR"
        # shellcheck disable=SC2086
        "${MMSEQS}" tar2db "${@}" "${TMP_PATH}/msa" --output-dbtype 11 --tar-include "${TAR2DB_INCLUDE}" ${THREADS_PAR} \
            || fail "tar2db died"
        rm -f -- "${@}"
        if [ -n "${SED_FIX_LOOKUP}" ]; then
            sed "${SED_FIX_LOOKUP}" "${TMP_PATH}/msa.lookup" > "${TMP_PATH}/msa.lookup_tmp"
            mv -f -- "${TMP_PATH}/msa.lookup_tmp" "${TMP_PATH}/msa.lookup"
        fi
        sed 's|\.raw_alg\.faa\.gz||g' "${TMP_PATH}/msa.lookup" > "${TMP_PATH}/msa.lookup.tmp"
        mv -f -- "${TMP_PATH}/msa.lookup.tmp" "${TMP_PATH}/msa.lookup"
        # shellcheck disable=SC2086
        "${MMSEQS}" msa2profile "${TMP_PATH}/msa" "${OUTDB}" --match-mode 1 --match-ratio 0.5 ${THREADS_PAR} \
            || fail "msa2profile died"
        mv -f -- "${TMP_PATH}/msa.lookup" "${OUTDB}.lookup"
        mv -f -- "${TMP_PATH}/msa.source" "${OUTDB}.source"
        if [ -n "${REMOVE_TMP}" ]; then
            # shellcheck disable=SC2086
            "${MMSEQS}" rmdb "${TMP_PATH}/msa" ${VERB_PAR} \
                || fail "rmdb died"
        fi
    ;;
    "GTDB")
        # shellcheck disable=SC2086
        "${MMSEQS}" tar2db "${TMP_PATH}/gtdb.tar.gz" "${TMP_PATH}/tardb" --tar-include 'faa.gz$' ${THREADS_PAR} \
            || fail "tar2db died"
        sed 's|_protein\.faa\.gz||g' "${TMP_PATH}/tardb.lookup" > "${TMP_PATH}/tardb.lookup.tmp"
        mv -f -- "${TMP_PATH}/tardb.lookup.tmp" "${TMP_PATH}/tardb.lookup"
        # shellcheck disable=SC2086
        "${MMSEQS}" createdb "${TMP_PATH}/tardb" "${OUTDB}" ${COMP_PAR} \
            || fail "createdb died"
        if [ -n "${REMOVE_TMP}" ]; then
            # shellcheck disable=SC2086
            "${MMSEQS}" rmdb "${TMP_PATH}/tardb" ${VERB_PAR} \
                || fail "rmdb died"
        fi
    ;;
esac
fi

if [ -n "${TAXONOMY}" ] && notExists "${OUTDB}_mapping"; then
    case "${SELECTION}" in
      "SILVA")
        mkdir -p "${TMP_PATH}/taxonomy"
        # shellcheck disable=SC2016
        CMD='BEGIN {
              ids["root"] = 1;
              print "1\t|\t1\t|\tno rank\t|\t-\t|" > taxdir"/nodes.dmp";
              print "1\t|\troot\t|\t-\t|\tscientific name\t|" > taxdir"/names.dmp";
          }
          {
              n = split($1, a, ";");
              gsub("domain", "superkingdom", $3);
              ids[$1] = $2;
              gsub(/[^,;]*;$/, "", $1);
              pname = $1;
              if (n == 2) {
                pname = "root";
              }
              pid = ids[pname];
              printf("%s\t|\t%s\t|\t%s\t|\t-\t|\n", $2, pid, $3) > taxdir"/nodes.dmp";
              printf("%s\t|\t%s\t|\t-\t|\tscientific name\t|\n", $2, a[n-1]) > taxdir"/names.dmp";
          }'
        awk -v taxdir="${TMP_PATH}/taxonomy" -F'\t' "$CMD" "${TMP_PATH}/silva_tax.txt"
        touch "${TMP_PATH}/taxonomy/merged.dmp"
        touch "${TMP_PATH}/taxonomy/delnodes.dmp"
        # shellcheck disable=SC2086
        "${MMSEQS}" createtaxdb "${OUTDB}" "${TMP_PATH}/taxdb" --ncbi-tax-dump "${TMP_PATH}/taxonomy" --tax-mapping-file "${TMP_PATH}/silva.acc_taxid" ${THREADS_PAR}
       ;;
     "NR")
        touch "${OUTDB}_mapping"
        # shellcheck disable=SC2086
        "${MMSEQS}" createtaxdb "${OUTDB}" "${TMP_PATH}/taxonomy" ${THREADS_PAR} \
            || fail "createtaxdb died"
        # shellcheck disable=SC2086
        "${MMSEQS}" nrtotaxmapping "${TMP_PATH}/pdb.accession2taxid" "${TMP_PATH}/prot.accession2taxid" "${OUTDB}" "${OUTDB}_mapping" ${THREADS_PAR} \
            || fail "nrtotaxmapping died"
       ;;
     "GTDB")
          # shellcheck disable=SC2016
          CMD='BEGIN {
              FS = "[\t;]";
              rank["c"] = "class";
              rank["d"] = "superkingdom";
              rank["f"] = "family";
              rank["g"] = "genus";
              rank["o"] = "order";
              rank["p"] = "phylum";
              rank["s"] = "species";
              taxCnt = 1;
              ids["root"] = 1;
              print "1\t|\t1\t|\tno rank\t|\t-\t|" > taxdir"/nodes.dmp";
              print "1\t|\troot\t|\t-\t|\tscientific name\t|" > taxdir"/names.dmp";
          }
          {
              prevTaxon=1;
              for (i = 2; i <= NF; i++) {
                  if ($i in ids) {
                      prevTaxon = ids[$i];
                  } else {
                      taxCnt++;
                      ids[$i] = taxCnt;
                      r = substr($i, 0, 1);
                      name = substr($i, 4);
                      gsub(/_/, " ", name);
                      printf("%s\t|\t%s\t|\t%s\t|\t-\t|\n", taxCnt, prevTaxon, rank[r]) > taxdir"/nodes.dmp";
                      printf("%s\t|\t%s\t|\t-\t|\tscientific name\t|\n", taxCnt, name) > taxdir"/names.dmp";
                      prevTaxon = taxCnt;
                  }
              }
              printf("%s\t%s\n", $1, ids[$NF]) > taxdir"/mapping_genomes";
          }'
          mkdir -p "${TMP_PATH}/taxonomy"
          awk -v taxdir="${TMP_PATH}/taxonomy" "$CMD" "${TMP_PATH}/bac120_taxonomy.tsv" "${TMP_PATH}/ar53_taxonomy.tsv"
          touch "${TMP_PATH}/taxonomy/merged.dmp"
          touch "${TMP_PATH}/taxonomy/delnodes.dmp"
          # shellcheck disable=SC2086
          "${MMSEQS}" createtaxdb "${OUTDB}" "${TMP_PATH}/taxdb" --ncbi-tax-dump "${TMP_PATH}/taxonomy" --tax-mapping-file "${TMP_PATH}/taxonomy/mapping_genomes" --tax-mapping-mode 1 ${THREADS_PAR} \
              || fail "createtaxdb died"
          if [ -n "${REMOVE_TMP}" ]; then
              rm -f -- "${TMP_PATH}/taxonomy/nodes.dmp" "${TMP_PATH}/taxonomy/names.dmp" "${TMP_PATH}/taxonomy/merged.dmp" "${TMP_PATH}/taxonomy/delnodes.dmp" "${TMP_PATH}/taxonomy/mapping_genomes" "${TMP_PATH}/bac120_taxonomy.tsv" "${TMP_PATH}/ar122_taxonomy.tsv"
              rm -rf -- "${TMP_PATH}/taxdb" "${TMP_PATH}/taxonomy"
          fi
     ;;
     *)
       # shellcheck disable=SC2086
       "${MMSEQS}" prefixid "${OUTDB}_h" "${TMP_PATH}/header_pref.tsv" --tsv ${THREADS_PAR} \
           || fail "prefixid died"
       awk '{ match($0, / OX=[0-9]+ /); if (RLENGTH != -1) { print $1"\t"substr($0, RSTART+4, RLENGTH-5); next; } match($0, / TaxID=[0-9]+ /); print $1"\t"substr($0, RSTART+7, RLENGTH-8); }' "${TMP_PATH}/header_pref.tsv" \
           | LC_ALL=C sort -n > "${OUTDB}_mapping"
       rm -f "${TMP_PATH}/header_pref.tsv"
       # shellcheck disable=SC2086
       "${MMSEQS}" createtaxdb "${OUTDB}" "${TMP_PATH}/taxonomy" ${THREADS_PAR} \
           || fail "createtaxdb died"
       ;;
     esac
fi

if notExists "${OUTDB}.version"; then
    mv -f "${TMP_PATH}/version" "${OUTDB}.version"
fi

if [ -n "${REMOVE_TMP}" ]; then
    rm -f "${TMP_PATH}/download.sh"
fi
