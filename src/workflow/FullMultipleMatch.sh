#!/bin/sh -e
# Aggregation Workflow
fail() {
    echo "Error: $1" | tee "logs.txt"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check amount of input variables
[ "$#" -ne 5 ] && echo "Please provide <QueryGenome> <TargetGenome> <Threads> <e-value Threshold> <p-value threshold>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;

RAWQUERYGENOME="$1"
RAWTARGETGENOME="$2"
THREADS="$3"
EVALUE="$4"
PVALUE="$5"
[ -f "logs.txt" ] && rm "logs.txt" ;

touch "logs.txt" ;

echo "Process started the : $(date)" >> "logs.txt"

####################### CREATEDB ############################

[ ! -d "Query_DB" ] &&  echo "Creation of Query_DB Directory" && mkdir -p "Query_DB";

[ ! -d "Target_DB" ] &&  echo "Creation of Target_DB Directory" && mkdir -p "Target_DB";


[ ! -d "Query_DB/Raw_Query_DB" ] && echo "Creation of Raw_Query_DB Directory" && mkdir -p "Query_DB/Raw_Query_DB" ;

if notExists "Query_DB/Raw_Query_DB/Query_Genome_DB" ; then
    echo "Creating DBs"
    echo "createdb of Query-set started at : $(date)" >> "logs.txt"

    /usr/bin/time -f "\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} createdb ${RAWQUERYGENOME} "Query_DB/Raw_Query_DB/Query_Genome_DB" \
    --dont-split-seq-by-len | tee -a "logs.txt" \
    || fail "createdb of query failed"
fi

[ ! -d "Target_DB/Raw_Target_DB" ] && echo "Creation of Raw_Target_DB Directory" && mkdir -p "Target_DB/Raw_Target_DB" ;



if notExists "Target_DB/Raw_Target_DB/Target_Genome_DB" ; then
    printf "\ncreatedb of Target-Set started at : $(date)" >> "logs.txt"

    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} createdb ${RAWTARGETGENOME} "Target_DB/Raw_Target_DB/Target_Genome_DB" \
     --dont-split-seq-by-len | tee -a "logs.txt" \
    || fail "createdb of target failed"
fi

########################### EXTRACTORFS ##########################

[ ! -d "Query_DB/Nucleotides_Orfs" ] && echo "Creation of Query_DB/Nucleotides_Orfs Directory" && mkdir -p "Query_DB/Nucleotides_Orfs";

if notExists "Query_DB/Nucleotides_Orfs/Nucl_Orfs" ; then
    echo "Extracting Orfs"
    printf "\nextractorfs of Query-Set started at : $(date)" >> "logs.txt"

    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} extractorfs "Query_DB/Raw_Query_DB/Query_Genome_DB"\
     "Query_DB/Nucleotides_Orfs/Nucl_Orfs" --min-length 30 --threads ${THREADS}| tee -a "logs.txt" \
    || fail "extractorfs on Query DB failed"
fi

[ ! -d "Target_DB/Nucleotides_Orfs" ] && echo "Creation of Target_DB/Nucleotides_Orfs Directory" && mkdir -p "Target_DB/Nucleotides_Orfs";

if notExists "Target_DB/Nucleotides_Orfs/Nucl_Orfs" ; then
    printf "\nextractorfs of Target-Set started at : $(date)" >> "logs.txt"

    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} extractorfs "Target_DB/Raw_Target_DB/Target_Genome_DB"\
     "Target_DB/Nucleotides_Orfs/Nucl_Orfs" --min-length 30 --threads ${THREADS}| tee -a "logs.txt" \
    || fail "extractorfs on Target DB failed"
fi

######################### TRANSLATENUCS ##################################

[ ! -d "Query_DB/AA_Orfs" ] && echo "Creation of Query_DB/AA_Orfs Directory" && mkdir -p "Query_DB/AA_Orfs" ;

if notExists "Query_DB/AA_Orfs/Orfs_AA" ; then
    echo "Translating Orfs"
    printf "\ntranslatenucs of Query Orfs started at : $(date)" >> "logs.txt"

    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} translatenucs "Query_DB/Nucleotides_Orfs/Nucl_Orfs" \
    "Query_DB/AA_Orfs/Orfs_AA" | tee -a "logs.txt" \
    || fail "translatenucs on Query"
fi

[ ! -d "Target_DB/AA_Orfs" ] && echo "Creation of Target_DB/AA_Orfs Directory" && mkdir -p "Target_DB/AA_Orfs" ;

if notExists "Target_DB/AA_Orfs/Orfs_AA" ; then
    printf "\ntranslatenucs of Target Orfs started at : $(date)" >> "logs.txt"

    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} translatenucs "Target_DB/Nucleotides_Orfs/Nucl_Orfs" \
    "Target_DB/AA_Orfs/Orfs_AA" | tee -a "logs.txt" \
    || fail "translatenucs on Target"
fi

####################################### NBR GENES ##########################################

if notExists "Query_DB/Nucleotides_Orfs/nbrORFsQuerySet" ; then

    # Only first arg is used
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} result2stats "Query_DB/Nucleotides_Orfs/Nucl_Orfs_set_lookup" \
    "Query_DB/Nucleotides_Orfs/Nucl_Orfs_set_lookup" "Query_DB/Nucleotides_Orfs/Nucl_Orfs_set_lookup" "Query_DB/Nucleotides_Orfs/nbrORFsQuerySet" \
    --stat linecount --threads ${THREADS} | tee -a "logs.txt" \
    || fail "result2stats for Query Set failed"
fi

if notExists "Target_DB/Nucleotides_Orfs/nbrORFsTargetSet" ; then
    # Only first arg is used
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} result2stats "Target_DB/Nucleotides_Orfs/Nucl_Orfs_set_lookup" \
    "Target_DB/Nucleotides_Orfs/Nucl_Orfs_set_lookup" "Target_DB/Nucleotides_Orfs/Nucl_Orfs_set_lookup" "Target_DB/Nucleotides_Orfs/nbrORFsTargetSet" \
    --stat linecount --threads ${THREADS} | tee -a "logs.txt" \
    || fail "result2stats for Target Set failed"
fi

####################################### RESEARCH #############################

[ ! -d "Main_Files" ] && echo "Creation of Main_Files Directory" && mkdir -p "Main_Files" ;

if notExists "Main_Files/1_SearchResults" ; then
    echo "Launching Research"
    printf "\nResearch Query VS Target started at : $(date)" >> "logs.txt"

    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} search "Query_DB/AA_Orfs/Orfs_AA" \
    "Target_DB/AA_Orfs/Orfs_AA" "Main_Files/1_SearchResults" "tmp" -e ${EVALUE} --threads ${THREADS} -s 7 -c 0.5 --max-seqs 1500 | tee -a "logs.txt" \
    || fail "search"
fi

####################################### ADD COLUMN ###########################################

if notExists "Main_Files/2_AddedColumn" ; then
    echo "Adding Column"
    printf "\nAdding column started at : $(date)" >> "logs.txt"
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} filterdb "Main_Files/1_SearchResults" \
    "Main_Files/2_AddedColumn" --join-db "Target_DB/Nucleotides_Orfs/Nucl_Orfs_orf_lookup" --column-to-take 0 | tee -a "logs.txt" \
    || fail "filterdb failed"
fi

############################# AGGREGATE BEST HIT ###############################
# REMOVE --simple-best-hit
if notExists "Main_Files/3_BestHitAggregation" ; then
    echo "BestHist Aggregation"
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} aggregate "Main_Files/2_AddedColumn"\
    "Main_Files/3_BestHitAggregation" "Target_DB/Nucleotides_Orfs/nbrORFsTargetSet" --set-column 10 --mode bestHit --threads ${THREADS} --simple-best-hit\
    | tee -a "logs.txt"\
    || fail "aggregate --mode bestHit failed"
fi

############################ MERGE BY QUERY SET ##################################

if notExists "Main_Files/4_Merged_BestHitAggregation" ; then
    echo "Merging BestHitAggreg"
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} mergeclusters "Main_Files/3_BestHitAggregation"\
    "Main_Files/4_Merged_BestHitAggregation" "Query_DB/Nucleotides_Orfs/Nucl_Orfs_set_lookup" --by-db --threads ${THREADS} | tee -a "logs.txt" \
    || fail "mergeclusters failed"
fi

##################### ADD HIT POSITION ON GENOME  #################################

if notExists "Main_Files/4.5_PositionsAdded" ; then
    echo "Adding Positions"
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} filterdb "Main_Files/4_Merged_BestHitAggregation" \
    "Main_Files/4.5_PositionsAdded" --compute-positions "Target_DB/Nucleotides_Orfs/Nucl_Orfs_orf_lookup" | tee -a "logs.txt" \
    || fail "Adding Postitions failed"
fi


######################### AGGREGATE BY PVAL ##############################

if notExists "Main_Files/5_PvalAggreg" ; then
    echo "Aggregate Pvals"
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} aggregate "Main_Files/4_Merged_BestHitAggregation" \
    "Main_Files/5_PvalAggreg" "Query_DB/Nucleotides_Orfs/nbrORFsQuerySet" "Target_DB/Nucleotides_Orfs/nbrORFsTargetSet" --mode pval \
    --threads ${THREADS} --set-column 10 --alpha ${PVALUE}| tee -a "logs.txt" \
    || fail "aggregate --mode pval failed"
fi

###################### GETTING HIT DISTANCES ###########################

if notExists "Main_Files/5.5_HitDistances" ; then
    echo "Assessing Hits distances"
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} aggregate "Main_Files/4.5_PositionsAdded" \
    "Main_Files/5.5_HitDistances" "Target_DB/Nucleotides_Orfs/nbrORFsTargetSet" "Query_DB/Nucleotides_Orfs/nbrORFsQuerySet" \
    "Target_DB/Raw_Target_DB/Target_Genome_DB" --mode clustering-index --threads ${THREADS} --alpha ${PVALUE} \
     --set-column 12 | tee -a "logs.txt" \
    || fail "Getting Hit Distances failed"
fi

###################### FORMATTING FILES TO TSV ##################################

if notExists "Main_Files/6_HeaderAddedDistances.tsv" ; then
    echo "Adding Header"
    /usr/bin/time -f "\t\t\t%E real,\t%U user,\t%S sys" -a -o "logs.txt" ${MMSEQS} createtsv "Query_DB/Raw_Query_DB/Query_Genome_DB" \
    "Target_DB/Raw_Target_DB/Target_Genome_DB" "Main_Files/5.5_HitDistances" "Main_Files/6_HeaderAddedDistances.tsv" --full-header\
    || fail "createtsv failed"
fi

if notExists "Main_Files/FinalData.tsv" ; then
    echo "Formating final File"

${MMSEQS} createtsv "Query_DB/Raw_Query_DB/Query_Genome_DB" \
    "Target_DB/Raw_Target_DB/Target_Genome_DB" "Main_Files/5_PvalAggreg" "Main_Files/5_PvalAggreg.tsv" --full-header

    cat "Main_Files/5_PvalAggreg" | tr -d "\0" | cut -f2 > "Main_Files/ToDelete" || fail "ToDelete"
    cat "Main_Files/6_HeaderAddedDistances.tsv" | tr -d "\0" > "Temp"
    paste "Temp" "Main_Files/ToDelete" > "Main_Files/FinalData.tsv" || fail "paste failed"
    rm "Main_Files/ToDelete" || fail "rm failed"
    rm "Temp"
fi

if notExists "Buoy" ; then
    echo "Creating Buoy"
    touch "Buoy"
fi

echo "Process ended the : $(date)" >> "logs.txt"

















