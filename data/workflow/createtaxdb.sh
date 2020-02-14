#!/bin/sh -e

notExists() {
	  [ ! -f "$1" ]
}

hasCommand() {
	  command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

notExists "$1" && echo "$1 not found!" && exit 1;

hasCommand wget
hasCommand awk
hasCommand zcat
hasCommand touch
hasCommand tar

TAXDBNAME="$1"
TMP_PATH="$2"

if [ "$DOWNLOAD_NCBITAXDUMP" -eq "1" ]; then
    # Download NCBI taxon information
    if notExists "${TMP_PATH}/ncbi_download.complete"; then
        echo "Download taxdump.tar.gz"
        wget -nv -O - "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" \
           | tar -C "${TMP_PATH}" -xzf - names.dmp nodes.dmp merged.dmp delnodes.dmp
        touch "${TMP_PATH}/ncbi_download.complete"
    fi
    NCBITAXINFO="${TMP_PATH}"
fi
if notExists "${TAXDBNAME}_mapping"; then
    if [ "$DOWNLOAD_MAPPING" -eq "1" ]; then
        # Download the latest UniProt ID mapping to extract taxon identifiers
        if notExists "${TMP_PATH}/mapping_download.complete"; then
            echo "Download idmapping.dat.gz"
            URL="https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
            wget -nv -O - "$URL" | zcat | awk '$2 == "NCBI_TaxID" {print $1"\t"$3 }' > "${TMP_PATH}/taxidmapping"
            touch "${TMP_PATH}/mapping_download.complete"
        fi
        MAPPINGFILE="${TMP_PATH}/taxidmapping"
    fi
    awk 'NR == FNR { f[$1] = $2; next } $2 in f { print $1"\t"f[$2] }' \
        "$MAPPINGFILE" "${TAXDBNAME}.lookup" > "${TAXDBNAME}_mapping"
fi

# finalize database
cp -f "${NCBITAXINFO}/names.dmp"     "${TAXDBNAME}_names.dmp"
cp -f "${NCBITAXINFO}/nodes.dmp"     "${TAXDBNAME}_nodes.dmp"
cp -f "${NCBITAXINFO}/merged.dmp"    "${TAXDBNAME}_merged.dmp"
cp -f "${NCBITAXINFO}/delnodes.dmp"  "${TAXDBNAME}_delnodes.dmp"
echo "Database created"

if [ -n "$REMOVE_TMP" ]; then
   rm -f "${TMP_PATH}/names.dmp" "${TMP_PATH}/nodes.dmp" "${TMP_PATH}/merged.dmp" "${TMP_PATH}/delnodes.dmp"
   rm -f "${TMP_PATH}/taxidmapping"
   if [ "$DOWNLOAD_DATA" -eq "1" ]; then
      rm -f "${TMP_PATH}/ncbi_download.complete" "${TMP_PATH}/mapping_download.complete"
   fi
   rm -f "${TMP_PATH}/createtaxdb.sh"
fi
