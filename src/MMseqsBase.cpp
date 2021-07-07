#include "Command.h"
#include "Parameters.h"
#include "CommandDeclarations.h"

Parameters& par = Parameters::getInstance();
std::vector<Command> baseCommands = {
        {"easy-search",          easysearch,           &par.easysearchworkflow,   COMMAND_EASY,
                "Sensitive homology search",
                "# Search multiple FASTA against FASTA (like BLASTP, TBLASTN, BLASTX, BLASTN --search-type 3, TBLASTX --search-type 2)\n"
                "mmseqs easy-search examples/QUERY.fasta examples/QUERY.fasta examples/DB.fasta result.m8 tmp\n\n"
                "# Iterative profile search from stdin (like PSI-BLAST)\n"
                "cat examples/QUERY.fasta | mmseqs easy-search stdin examples/DB.fasta result.m8 tmp --num-iterations 2\n\n"
                "# Profile search against small databases (e.g. PFAM, eggNOG)\n"
                "mmseqs databases PFAM pfam_db tmp\n"
                "mmseqs easy-search examples/QUERY.fasta pfam_db res.m8 tmp\n\n"
                "# Exhaustive search against sequences or profiles (works for large DBs)\n"
                "mmseqs easy-search examples/QUERY.fasta targetProfiles res.m8 tmp --exhaustive-search\n\n"
                "# Increasing sensitivity search (from 2 to 7 in 3 steps)\n"
                "mmseqs easy-search examples/QUERY.fasta examples/DB.fasta result.m8 tmp --start-sens 2 -s 7 --sens-steps 3\n",
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryFastaFile1[.gz|.bz2]> ... <i:queryFastaFileN[.gz|.bz2]>|<i:stdin> <i:targetFastaFile[.gz]>|<i:targetDB> <o:alignmentFile> <tmpDir>",
                CITATION_SERVER | CITATION_MMSEQS2,{{"fastaFile[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &DbValidator::flatfileAndStdin },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-linsearch",       easylinsearch,        &par.easylinsearchworkflow,COMMAND_EASY | COMMAND_EXPERT,
                "Fast, less sensitive homology search",
                NULL,
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryFastaFile1[.gz|.bz2]> ... <i:queryFastaFileN[.gz|.bz2]> <i:targetFastaFile[.gz|.bz2]>|<i:targetDB> <o:alignmentFile> <tmpDir>",
                CITATION_MMSEQS2|CITATION_LINCLUST, {{"fastaFile[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &DbValidator::flatfileAndStdin },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-cluster",         easycluster,          &par.easyclusterworkflow, COMMAND_EASY,
                "Slower, sensitive clustering",
                "mmseqs easy-cluster examples/DB.fasta result tmp\n"
                "# Cluster output\n"
                "#  - result_rep_seq.fasta: Representatives\n"
                "#  - result_all_seq.fasta: FASTA-like per cluster\n"
                "#  - result_cluster.tsv:   Adjacency list\n\n"
                "# Important parameter: --min-seq-id, --cov-mode and -c \n"
                "#                  --cov-mode \n"
                "#                  0    1    2\n"
                "# Q: MAVGTACRPA  60%  IGN  60%\n"
                "# T: -AVGTAC---  60% 100%  IGN\n"
                "#        -c 0.7    -    +    -\n"
                "#        -c 0.6    +    +    +\n\n"
                "# Cascaded clustering with reassignment\n"
                "# - Corrects criteria-violoations of cascaded merging\n"
                "# - Produces more clusters and is a bit slower\n"
                "mmseqs easy-cluster examples/DB.fasta result tmp --cluster-reassign\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:fastaFile1[.gz|.bz2]> ... <i:fastaFileN[.gz|.bz2]> <o:clusterPrefix> <tmpDir>",
                CITATION_MMSEQS2|CITATION_LINCLUST, {{"queryFastaFile[.gz]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &DbValidator::flatfileAndStdin },
                                                           {"clusterPrefix", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-linclust",        easylinclust,         &par.easylinclustworkflow, COMMAND_EASY,
                "Fast linear time cluster, less sensitive clustering",
                "mmseqs easy-linclust examples/DB.fasta result tmp\n\n"
                "# Linclust output\n"
                "#  - result_rep_seq.fasta: Representatives\n"
                "#  - result_all_seq.fasta: FASTA-like per cluster\n"
                "#  - result_cluster.tsv:   Adjecency list\n\n"
                "# Important parameter: --min-seq-id, --cov-mode and -c \n"
                "#                  --cov-mode \n"
                "#                  0    1    2\n"
                "# Q: MAVGTACRPA  60%  IGN  60%\n"
                "# T: -AVGTAC---  60% 100%  IGN\n"
                "#        -c 0.7    -    +    -\n"
                "#        -c 0.6    +    +    +\n\n"
                "# Cluster nucleotide sequences \n"
                "mmseqs easy-linclust examples/DB.fasta result tmp --kmer-per-seq-scale 0.3\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:fastaFile1[.gz|.bz2]> ... <i:fastaFileN[.gz|.bz2]> <o:clusterPrefix> <tmpDir>",
                CITATION_MMSEQS2|CITATION_LINCLUST, {{"fastaFile[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &DbValidator::flatfileAndStdin },
                                                            {"clusterPrefix", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                            {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-taxonomy",        easytaxonomy,         &par.easytaxonomy,         COMMAND_EASY,
                "Taxonomic classification",
                "# Assign taxonomic labels to FASTA sequences\n"
                "  - result_tophit_aln: top hits\n"
                "  - result_tophit_report: coverage profiles per database entry\n"
                "  - result_report: kraken style report\n"
                "# Download a sequence database with taxonomy information\n"
                "mmseqs databases UniProtKB/Swiss-Prot swissprotDB tmp\n\n"
                "# Assign taxonomy based on 2bLCA hit\n"
                "mmseqs easy-taxonomy examples/DB.fasta swissprotDB result tmp\n\n"
                "# Assign taxonomy based on top hit\n"
                "mmseqs easy-taxonomy examples/DB.fasta swissprotDB result tmp --lca-mode 4\n\n"
                "# Assign taxonomy without ORF prefilter\n"
                "# Classifies higher percentage for short nucleotide input (e.g. short reads) at the cost of speed\n"
                "mmseqs easy-taxonomy queryNuclDB swissprotDB result tmp --orf-filter 0\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:fastaFile1[.gz|.bz2]> ... <i:fastaFileN[.gz|.bz2]> <i:targetDB> <o:taxReports> <tmpDir>",
                CITATION_TAXONOMY|CITATION_MMSEQS2, {{"queryFastaFile[.gz]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_DATA|DbType::VARIADIC, &DbValidator::flatfileAndStdin },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                                           {"taxReports",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"easy-rbh",                  easyrbh,                  &par.easysearchworkflow,       COMMAND_EASY,
                "Find reciprocal best hit",
                "# Assign reciprocal best hit\n"
                "mmseqs easy-rbh examples/QUERY.fasta examples/DB.fasta result tmp\n\n",
                "Eli Levy Karin & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryFastaFile1[.gz|.bz2]> <i:targetFastaFile[.gz|.bz2]>|<i:targetDB> <o:alignmentFile> <tmpDir>",
                CITATION_MMSEQS2,{{"fastaFile[.gz|.bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &DbValidator::flatfileAndStdin },
                                   {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"databases",            databases,            &par.databases,            COMMAND_DATABASE_CREATION,
                "List and download databases",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<name> <o:sequenceDB> <tmpDir>",
                CITATION_TAXONOMY|CITATION_MMSEQS2, {{"selection", 0, DbType::ZERO_OR_ALL, &DbValidator::empty },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"tmpDir",     DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"createdb",             createdb,             &par.createdb,             COMMAND_DATABASE_CREATION,
                "Convert FASTA/Q file(s) to a sequence DB",
                "# Create a sequence database from multiple FASTA files\n"
                "mmseqs createdb file1.fa file2.fa.gz file3.fa sequenceDB\n\n"
                "# Create a seqDB from stdin\n"
                "cat seq.fasta | mmseqs createdb stdin sequenceDB\n\n"
                "# Create a seqDB by indexing existing FASTA/Q (for single line fasta entries only)\n"
                "mmseqs createdb seq.fasta sequenceDB --createdb-mode 1\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:fastaFile1[.gz|.bz2]> ... <i:fastaFileN[.gz|.bz2]>|<i:stdin> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"fast[a|q]File[.gz|bz2]|stdin", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfileStdinAndGeneric },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"indexdb",              indexdb,              &par.indexdb,              COMMAND_HIDDEN,
                NULL,
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <o:sequenceIndexDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                                           {"sequenceIndexDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"createindex",          createindex,          &par.createindex,          COMMAND_DATABASE_CREATION,
                "Store precomputed index on disk to reduce search overhead",
                "# Create protein sequence index\n"
                "mmseqs createindex sequenceDB tmp\n\n"
                "# Create TBLASTX/N index from nucleotide sequences\n"
                "mmseqs createindex sequenceDB tmp --search-type 2\n\n"
                "# Create BLASTN index from nucleotide sequences\n"
                "mmseqs createindex sequenceDB tmp --search-type 3\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <tmpDir>",
                CITATION_SERVER | CITATION_MMSEQS2,{{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"createlinindex",       createlinindex,       &par.createlinindex,       COMMAND_DATABASE_CREATION | COMMAND_EXPERT,
                "Create linsearch index",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <tmpDir>",
                CITATION_SERVER | CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"convertmsa",           convertmsa,           &par.convertmsa,           COMMAND_DATABASE_CREATION,
                "Convert Stockholm/PFAM MSA file to a MSA DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:msaFile.sto[.gz]> <o:msaDB>",
                CITATION_SERVER |CITATION_MMSEQS2, {{"msaFile.sto[.gz]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"msaDB",DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::msaDb }}},
        {"tsv2db",               tsv2db,               &par.tsv2db,               COMMAND_DATABASE_CREATION | COMMAND_EXPERT,
                "Convert a TSV file to any DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:tsvFile> <o:resultDB>",
                CITATION_MMSEQS2, {{"tsvFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                          {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
        {"tar2db",               tar2db,               &par.tar2db,               COMMAND_DATABASE_CREATION | COMMAND_EXPERT,
                "Convert content of tar archives to any DB",
                "# Assuming tar archive containing three aligned FASTA files:\n"
                "#  * folder/msa1.fa.gz  * folder/msa2.fa  * folder/msa3.fa\n"
                "# Create a msaDB with three DB entries each containing a separate MSA\n"
                "mmseqs tar2db archive.tar.gz msaDB --output-dbtype 11\n",
                "Milot Mirdita <milot@mirdita.de>",
                "<i:tar[.gz]> ... <i:tar[.gz]> <o:resultDB>",
                CITATION_MMSEQS2, {{".tar[.gz]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},


        {"search",               search,               &par.searchworkflow,       COMMAND_MAIN,
                "Sensitive homology search",
                "# Search multiple FASTA against FASTA (like BLASTP, TBLASTN, BLASTX, BLASTN --search-type 3, TBLASTX --search-type 2)\n"
                "mmseqs search queryDB targetDB resultDB tmp\n"
                "mmseqs convertalis queryDB targetDB resultDB result.m8\n\n"
                "# Iterative profile search (like PSI-BLAST)\n"
                "mmseqs search queryDB targetDB resultDB tmp --num-iterations 2\n\n"
                "# Profile search against small databases (e.g. PFAM, eggNOG)\n"
                "mmseqs databases PFAM pfam_db tmp\n"
                "mmseqs search queryDB pfam_db resultDB tmp\n\n"
                "# Exhaustive search against sequences or profiles (works for large DBs)\n"
                "mmseqs search queryDB targetDB resultDB tmp --exhaustive-search\n\n"
                "# Increasing sensitivity search (from 2 to 7 in 3 steps)\n"
                "mmseqs search queryDB targetDB resultDB --start-sens 2 -s 7 --sens-steps 3\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"linsearch",            linsearch,            &par.linsearchworkflow,    COMMAND_MAIN|COMMAND_EXPERT,
                "Fast, less sensitive homology search",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"map",                  map,                  &par.mapworkflow,          COMMAND_MAIN,
                "Map nearly identical sequences",
                NULL,
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_PLASS|CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"rbh",                  rbh,                  &par.searchworkflow,       COMMAND_MAIN,
                "Reciprocal best hit search",
                NULL,
                "Eli Levy Karin",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"linclust",          linclust,          &par.linclustworkflow,           COMMAND_MAIN,
                "Fast, less sensitive clustering",
                "# Linear-time clustering of FASTA file\n"
                "mmseqs linclust sequenceDB clusterDB tmp\n\n"
                "                   --cov-mode \n"
                "# Sequence         0    1    2\n"
                "# Q: MAVGTACRPA  60%  IGN  60%\n"
                "# T: -AVGTAC---  60% 100%  IGN\n"
                "# Cutoff -c 0.7    -    +    -\n"
                "#        -c 0.6    +    +    +\n\n"
                "# Cluster nucleotide sequences \n"
                "mmseqs easy-linclust nucl.fasta result tmp --kmer-per-seq-scale 0.3\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <o:clusterDB> <tmpDir>",
                CITATION_MMSEQS2|CITATION_LINCLUST, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"cluster",              clusteringworkflow,   &par.clusterworkflow,      COMMAND_MAIN,
                "Slower, sensitive clustering",
                "# Cascaded clustering of FASTA file\n"
                "mmseqs cluster sequenceDB clusterDB tmp\n\n"
                "#                  --cov-mode \n"
                "# Sequence         0    1    2\n"
                "# Q: MAVGTACRPA  60%  IGN  60%\n"
                "# T: -AVGTAC---  60% 100%  IGN\n"
                "# Cutoff -c 0.7    -    +    -\n"
                "#        -c 0.6    +    +    +\n\n"
                "# Cascaded clustering with reassignment\n"
                "# - Corrects criteria-violoations of cascaded merging\n"
                "# - Produces more clusters and is a bit slower\n"
                "mmseqs cluster sequenceDB clusterDB tmp --cluster-reassign\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr> & Lars von den Driesch",
                "<i:sequenceDB> <o:clusterDB> <tmpDir>",
                CITATION_LINCLUST|CITATION_MMSEQS1|CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"clusterupdate",        clusterupdate,        &par.clusterUpdate,        COMMAND_MAIN,
                "Update previous clustering with new sequences",
                "# Update clustering workflow \n"
                "# Perform initial clustering of 5000 sequences\n"
                "mmseqs createdb <(head -n 10000 examples/DB.fasta) sequenceDB\n"
                "mmseqs cluster sequenceDB clusterDB tmp\n\n"
                "# Use-case 1: Update by only adding sequences\n"
                "mmseqs createdb examples/QUERY.fasta addedSequenceDB\n"
                "mmseqs concatdbs sequenceDB addedSequenceDB allSequenceDB\n"
                "mmseqs concatdbs sequenceDB_h addedSequenceDB_h allSequenceDB_h\n"
                "mmseqs clusterupdate sequenceDB allSequenceDB clusterDB newSequenceDB newClusterDB tmp\n\n"
                "# Use-case 2: Update clustering with deletions)\n"
                "# Create a FASTA file missing 500 of the original sequences and 2500 new ones\n"
                "mmseqs createdb <(tail -n +1001 examples/DB.fasta | head -n 15000) updateSequenceDB\n"
                "mmseqs clusterupdate sequenceDB updateSequenceDB clusterDB newSequenceDB newClusterDB tmp\n",
                "Clovis Galiez & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:oldSequenceDB> <i:newSequenceDB> <i:oldClustResultDB> <o:newMappedSequenceDB> <o:newClustResultDB> <tmpDir>",
                CITATION_MMSEQS2|CITATION_MMSEQS1,{{"oldSequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER|DbType::NEED_LOOKUP, &DbValidator::sequenceDb },
                                                          {"newSequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER|DbType::NEED_LOOKUP, &DbValidator::sequenceDb },
                                                          {"oldClustResultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                                                          {"newMappedSequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                                          {"newClustResultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb},
                                                          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"taxonomy",             taxonomy,             &par.taxonomy,             COMMAND_MAIN,
                "Taxonomic classification",
                "# Download a sequence database with taxonomy information\n"
                "mmseqs databases UniProtKB/Swiss-Prot swissprotDB tmp\n\n"
                "# Assign taxonomy based on 2bLCA\n"
                "mmseqs taxonomy queryDB swissprotDB result tmp\n\n"
                "# Assign taxonomy based on top hit\n"
                "mmseqs taxonomy queryDB swissprotDB result tmp --lca-mode 4\n\n"
                "# Assign taxonomy without ORF prefilter\n"
                "# Classifies higher percentage for short nucleotide input (e.g. short reads) at the cost of speed\n"
                "mmseqs taxonomy queryNuclDB swissprotDB result tmp --orf-filter 0\n\n"
                "# Create a Krona report\n"
                "mmseqs taxonomyreport swissprotDB result report.html --report-mode 1\n",
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr> & Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:queryDB> <i:targetDB> <o:taxaDB> <tmpDir>",
                CITATION_TAXONOMY|CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                                          {"taxaDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxResult },
                                                          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},



        {"convertalis",          convertalignments,    &par.convertalignments,    COMMAND_FORMAT_CONVERSION,
                "Convert alignment DB to BLAST-tab, SAM or custom format",
                "# Create output in BLAST M8 format (12 columns):\n"
                "#  (1,2) identifiers for query and target sequences/profiles,\n"
                "#  (3) sequence identity, (4) alignment length, (5) number of mismatches,\n"
                "#  (6) number of gap openings, (7-8, 9-10) alignment start and end-position in query and in target,\n"
                "#  (11) E-value, and (12) bit score\n"
                "mmseqs convertalis queryDB targetDB result.m8\n\n"
                "# Create a TSV containing pairwise alignments\n"
                "mmseqs convertalis queryDB targetDB result.tsv --format-output query,target,qaln,taln\n\n"
                "# Annotate a alignment result with taxonomy information from targetDB\n"
                "mmseqs convertalis queryDB targetDB result.tsv --format-output query,target,taxid,taxname,taxlineage\n\n"
                " Create SAM output\n"
                "mmseqs convertalis queryDB targetDB result.sam --format-mode 1\n\n"
                "# Create a TSV containing which query file a result comes from\n"
                "mmseqs createdb euk_queries.fasta bac_queries.fasta queryDB\n"
                "mmseqs convertalis queryDB targetDB result.tsv --format-output qset,query,target\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDb> <i:targetDb> <i:alignmentDB> <o:alignmentFile>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"alignmentFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"createtsv",            createtsv,            &par.createtsv,            COMMAND_FORMAT_CONVERSION,
                "Convert result DB to tab-separated flat file",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> [<i:targetDB>] <i:resultDB> <o:tsvFile>",
                CITATION_MMSEQS2,{{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"convert2fasta",        convert2fasta,        &par.convert2fasta,        COMMAND_FORMAT_CONVERSION,
                "Convert sequence DB to FASTA format",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:sequenceDB> <o:fastaFile>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::allDb },
                                                           {"fastaFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"result2flat",          result2flat,          &par.result2flat,          COMMAND_FORMAT_CONVERSION | COMMAND_EXPERT,
                "Create flat file by adding FASTA headers to DB entries",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:fastaDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                    {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_HEADER, &DbValidator::sequenceDb },
                    {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                    {"fastaDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"createseqfiledb",      createseqfiledb,      &par.createseqfiledb,      COMMAND_FORMAT_CONVERSION | COMMAND_EXPERT,
                "Create a DB of unaligned FASTA entries",
                "# Gather all sequences from a cluster DB\n"
                "mmseqs createseqfiledb sequenceDB clusterDB unalignedDB --min-sequences 2\n"
                "# Build MSAs with Clustal-Omega\n"
                "mmseqs apply unalignedDB msaDB -- clustalo -i - -o stdout --threads=1\n",
                "Milot Mirdita <milot@mirdita.de>",
                "<i:sequenceDB> <i:resultDB> <o:fastaDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"fastaDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb }}},



        {"createtaxdb",          createtaxdb,          &par.createtaxdb,          COMMAND_TAXONOMY,
                "Add taxonomic labels to sequence DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <tmpDir>",
                CITATION_TAXONOMY, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"createbintaxonomy",    createbintaxonomy,    &par.onlyverbosity,        COMMAND_TAXONOMY | COMMAND_EXPERT,
                "Create binary taxonomy from NCBI input",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:names.dmp> <i:nodes.dmp> <i:merged.dmp> <o:taxonomyFile>",
                CITATION_TAXONOMY, {{"names.dmp", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"nodes.dmp", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"merged.dmp", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"taxonomyFile",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"addtaxonomy",          addtaxonomy,          &par.addtaxonomy,          COMMAND_TAXONOMY | COMMAND_EXPERT,
                "Add taxonomic labels to result DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:targetDB> <i:resultDB> <o:resultDB>",
                CITATION_TAXONOMY, {{"targetDB", DbType::ACCESS_MODE_INPUT|DbType::NEED_TAXONOMY, DbType::NEED_DATA, &DbValidator::taxSequenceDb },
                                                           {"resultDB",   DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                                           {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"taxonomyreport",       taxonomyreport,       &par.taxonomyreport,       COMMAND_TAXONOMY | COMMAND_FORMAT_CONVERSION,
                "Create a taxonomy report in Kraken or Krona format",
                NULL,
                "Milot Mirdita <milot@mirdita.de> & Florian Breitwieser <florian.bw@gmail.com>",
                "<i:seqTaxDB> <i:taxResultDB/resultDB/sequenceDB> <o:taxonomyReport>",
                CITATION_TAXONOMY, {{"seqTaxDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                                           {"taxResultDB/resultDB/sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC,  &DbValidator::taxonomyReportInput },
                                                           {"taxonomyReport",    DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"filtertaxdb",          filtertaxdb,          &par.filtertaxdb,          COMMAND_TAXONOMY,
                "Filter taxonomy result database",
                "# Download a sequence database with taxonomy information\n"
                "mmseqs databases UniProtKB/Swiss-Prot swissprotDB tmp\n"
                "# Annotate a queryDB with taxonomy information\n"
                "mmseqs taxonomy queryDB swissprotDB taxDB tmp\n\n"
                "# Retain all unclassified hits\n"
                "mmseqs filtertaxdb swissprotDB taxDB filteredTaxDB --taxon-list 0\n"
                "mmseqs createsubdb <(awk '$3 == 1' filteredTaxDB.index) queryDB queryUnclassifiedDB\n\n"
                "# Retain all eukaryotic hits except fungi\n"
                "mmseqs filtertaxdb swissprotDB taxDB filteredTaxDB --taxon-list '2759&&!4751'\n\n"
                "# Retain all human and chlamydia hits\n"
                "mmseqs filtertaxdb swissprotDB taxDB filteredTaxDB --taxon-list '9606||810'\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:targetDB> <i:taxDB> <o:taxDB>",
                CITATION_TAXONOMY, {{"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::taxResult },
                                                           {"taxDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxResult }}},
        // TODO make consistent with seqTaxDB -> taxSeqDb in Wiki
        {"filtertaxseqdb",       filtertaxseqdb,       &par.filtertaxseqdb,       COMMAND_TAXONOMY,
                "Filter taxonomy sequence database",
                "# Download a sequence database with taxonomy information\n"
                "mmseqs databases UniProtKB/Swiss-Prot swissprotDB tmp\n\n"
                "# Retain all bacterial sequences\n"
                "mmseqs filtertaxseqdb swissprotDB swissprotDB_only_bac --taxon-list 2\n\n"
                "# Retain all eukaryotic sequences except fungi\n"
                "mmseqs filtertaxseqdb swissprotDB swissprotDB_euk_wo_fungi --taxon-list '2759&&!4751'\n\n"
                "# Retain all human and chlamydia sequences\n"
                "mmseqs filtertaxseqdb swissprotDB swissprotDB_human_and_chlamydia --taxon-list '9606||810'\n\n",
                "Eli Levy Karin <eli.levy.karin@gmail.com> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:taxSeqDB> <o:taxSeqDB>",
                CITATION_TAXONOMY, {{"taxSeqDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                                           {"taxSeqDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxSequenceDb }}},
        {"aggregatetax",         aggregatetax,         &par.aggregatetax,         COMMAND_TAXONOMY,
                "Aggregate multiple taxon labels to a single label",
                "# Download a sequence database with taxonomy information\n"
                "mmseqs databases UniProtKB/Swiss-Prot swissprotDB tmp\n\n"
                "# Create a nucleotide sequence database from FASTA\n"
                "mmseqs createdb contigs.fasta contigsDb\n\n"
                "# Extract all orfs from each contig and translate them\n"
                "mmseqs extractorfs contigsDb orfsAaDb --translate\n\n"
                "# Assign taxonomy to each orf\n"
                "mmseqs taxonomy orfsAaDb swissprotDB taxPerOrf tmp\n\n"
                "# Aggregate taxonomic assignments on each contig\n"
                "mmseqs aggregatetax swissprotDB orfsAaDb_h taxPerOrf taxPerContig --majority 0.5\n\n",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:taxSeqDB> <i:setToSeqMap> <i:taxResPerSeqDB> <o:taxResPerSetDB>",
                CITATION_TAXONOMY, {{"taxSeqDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                                           {"setToSeqMap",   DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                                           {"taxResPerSeqDB",   DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::taxResult },
                                                           {"taxResPerSetDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxResult }}},
        {"aggregatetaxweights",  aggregatetaxweights,  &par.aggregatetaxweights,  COMMAND_TAXONOMY,
                "Aggregate multiple taxon labels to a single label",
                "# Download a sequence database with taxonomy information\n"
                "mmseqs databases UniProtKB/Swiss-Prot swissprotDB tmp\n\n"
                "# Create a nucleotide sequence database from FASTA\n"
                "mmseqs createdb contigs.fasta contigsDb\n\n"
                "# Extract all orfs from each contig and translate them\n"
                "mmseqs extractorfs contigsDb orfsAaDb --translate\n\n"
                "# Assign taxonomy to each orf\n"
                "mmseqs taxonomy orfsAaDb swissprotDB taxPerOrf tmp --tax-output-mode 2\n\n"
                "# Aggregate taxonomic assignments on each contig\n"
                "mmseqs aggregatetaxweights swissprotDB orfsAaDb_h taxPerOrf taxPerOrf_aln taxPerContig --majority 0.5\n\n",
                "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:taxSeqDB> <i:setToSeqMap> <i:taxResPerSeqDB> <i:taxAlnResPerSeqDB> <o:taxResPerSetDB>",
                CITATION_TAXONOMY, {{"taxSeqDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                                           {"setToSeqMap",   DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                                           {"taxResPerSeqDB",   DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::taxResult },
                                                           {"taxAlnResPerSeqDB",   DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"taxResPerSetDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxResult }}},
        {"lcaalign",             lcaalign,             &par.align,                COMMAND_TAXONOMY,
                "Efficient gapped alignment for lca computation",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:alignmentDB>",
                CITATION_TAXONOMY, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::taxSequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"lca",                  lca,                  &par.lca,                  COMMAND_TAXONOMY,
                "Compute the lowest common ancestor",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:targetDB> <i:resultDB> <o:taxaDB>",
                CITATION_TAXONOMY|CITATION_MMSEQS2, {{"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"taxDB",    DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxResult }}},
        {"majoritylca",          majoritylca,          &par.majoritylca,          COMMAND_TAXONOMY | COMMAND_EXPERT,
                "Compute the lowest common ancestor using majority voting",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:targetDB> <i:resultDB> <o:taxaDB>",
                CITATION_TAXONOMY, {{"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"taxDB",    DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxResult }}},


        {"multihitdb",           multihitdb,           &par.multihitdb,           COMMAND_MULTIHIT,
                "Create sequence DB for multi hit searches",
                NULL,
                "Ruoshi Zhang, Clovis Norroy & Milot Mirdita <milot@mirdita.de>",
                "<i:fastaFile1[.gz|bz2]> ... <i:fastaFileN[.gz|bz2]> <o:setDB> <tmpDir>",
                CITATION_MMSEQS2, {{"fast[a|q]File[.gz|bz2]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                                           {"setDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"multihitsearch",       multihitsearch,       &par.multihitsearch,       COMMAND_MULTIHIT,
                "Search with a grouped set of sequences against another grouped set",
                NULL,
                "Ruoshi Zhang, Clovis Norroy & Milot Mirdita <milot@mirdita.de>",
                "<i:querySetDB> <i:targetSetDB> <o:resultDB> <tmpDir>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"besthitperset",        besthitperset,        &par.besthitbyset,         COMMAND_MULTIHIT,
                "For each set of sequences compute the best element and update p-value",
                NULL,
                "Ruoshi Zhang, Clovis Norroy & Milot Mirdita <milot@mirdita.de>",
                " <i:targetSetDB> <i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
        {"combinepvalperset",    combinepvalperset,    &par.combinepvalbyset,     COMMAND_MULTIHIT,
                "For each set compute the combined p-value",
                NULL,
                "Ruoshi Zhang, Clovis Norroy & Milot Mirdita <milot@mirdita.de>",
                "<i:querySetDB> <i:targetSetDB> <i:resultDB> <o:pvalDB>",
                CITATION_MMSEQS2, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"pvalDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"mergeresultsbyset",    mergeresultsbyset,    &par.threadsandcompression,COMMAND_MULTIHIT,
                "Merge results from multiple ORFs back to their respective contig",
                NULL,
                "Ruoshi Zhang, Clovis Norroy & Milot Mirdita <milot@mirdita.de>",
                "<i:setDB> <i:DB> <o:DB>",
                CITATION_MMSEQS2, {{"setDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                                           {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},



        {"prefilter",            prefilter,            &par.prefilter,            COMMAND_PREFILTER,
                "Double consecutive diagonal k-mer search",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> & Maria Hauser",
                "<i:queryDB> <i:targetDB> <o:prefilterDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"prefilterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::prefilterDb }}},

        {"ungappedprefilter",    ungappedprefilter,    &par.ungappedprefilter,    COMMAND_PREFILTER,
                "Optimal diagonal score search",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <o:prefilterDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"prefilterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::prefilterDb }}},
        {"kmermatcher",          kmermatcher,          &par.kmermatcher,          COMMAND_PREFILTER,
                "Find bottom-m-hashed k-mer matches within sequence DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <o:prefilterDB>",
                CITATION_MMSEQS2,{{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"prefilterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::prefilterDb }}},
        {"kmersearch",           kmersearch,           &par.kmersearch,           COMMAND_PREFILTER,
                "Find bottom-m-hashed k-mer matches between target and query DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:kmerIndexDB> <o:prefilterDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::indexDb },
                                         {"prefilterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::prefilterDb }}},
        {"kmerindexdb",          kmerindexdb,          &par.kmerindexdb,          COMMAND_HIDDEN,
                "Create bottom-m-hashed k-mer index",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> ",
                "<i:sequenceDB> <o:kmerIndexDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"kmerIndexDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},



        {"align",                align,                &par.align,                COMMAND_ALIGNMENT,
                "Optimal gapped local alignment",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> & Maria Hauser",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:alignmentDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"alignall",             alignall,             &par.alignall,             COMMAND_ALIGNMENT,
                "Within-result all-vs-all gapped local alignment",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <i:resultDB> <o:alignmentDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"transitivealign",      transitivealign,      &par.align,                COMMAND_ALIGNMENT,
                "Transfer alignments via transitivity",
                //"Infer the alignment A->C via B, B being the center sequence and A,C each pairwise aligned against B",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <i:alignmentDB> <o:alignmentDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                                           {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"rescorediagonal",     rescorediagonal,       &par.rescorediagonal,      COMMAND_ALIGNMENT,
                "Compute sequence identity for diagonal",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:prefilterDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"alignbykmer",         alignbykmer,           &par.alignbykmer,          COMMAND_ALIGNMENT,
                "Heuristic gapped local k-mer based alignment",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},



        {"clust",                clust,                &par.clust,                COMMAND_CLUSTER,
                "Cluster result by Set-Cover/Connected-Component/Greedy-Incremental",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> & Lars von den Driesch & Maria Hauser",
                "<i:sequenceDB> <i:resultDB> <o:clusterDB>",
                CITATION_MMSEQS2|CITATION_MMSEQS1,{{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                          {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb }}},
        {"clusthash",            clusthash,            &par.clusthash,            COMMAND_CLUSTER,
                "Hash-based clustering of equal length sequences",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> ",
                "<i:sequenceDB> <o:alignmentDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"mergeclusters",        mergeclusters,        &par.threadsandcompression,COMMAND_CLUSTER,
                "Merge multiple cascaded clustering steps",
                NULL,
                "Maria Hauser & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <o:clusterDB> <i:clusterDB1> ... <i:clusterDBn>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                          {"clusterDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                                                          {"clusterDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::clusterDb }}},



        {"compress",             compress,             &par.onlythreads,          COMMAND_STORAGE,
                "Compress DB entries",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:DB> <o:DB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                                           {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"decompress",           decompress,           &par.onlythreads,          COMMAND_STORAGE,
                "Decompress DB entries",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:DB> <o:DB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                                           {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"rmdb",                 rmdb,                 &par.onlyverbosity,        COMMAND_STORAGE,
                "Remove a DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:DB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL }}},
        {"mvdb",                 mvdb,                 &par.onlyverbosity,        COMMAND_STORAGE,
                "Move a DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:srcDB> <o:dstDB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"cpdb",                 cpdb,                 &par.onlyverbosity,        COMMAND_STORAGE,
                "Copy a DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:srcDB> <o:dstDB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"lndb",                 lndb,                 &par.onlyverbosity,        COMMAND_STORAGE,
                "Symlink a DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:srcDB> <o:dstDB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"unpackdb",             unpackdb,             &par.onlyverbosity,        COMMAND_STORAGE,
                "Unpack a DB into separate files",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:DB> <o:outDir>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL },
                                          {"outDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"touchdb",              touchdb,              &par.onlythreads,          COMMAND_STORAGE,
                "Preload DB into memory (page cache)",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> ",
                "<i:DB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb }}},


        {"createsubdb",          createsubdb,          &par.createsubdb,          COMMAND_SET,
                "Create a subset of a DB from list of DB keys",
                "# Create a new sequenceDB from sequenceDB entries with keys 1, 2 and 3\n"
                "mmseqs createsubdb <(printf '1\n2\n3\n') sequenceDB oneTwoThreeDB\n\n"
                "# Create a new sequence database with representatives of clusterDB\n"
                "mmseqs cluster sequenceDB clusterDB tmp\n"
                "mmseqs createsubdb clusterDB sequenceDB representativesDB\n",
                "Milot Mirdita <milot@mirdita.de>",
                "<i:subsetFile|DB> <i:DB> <o:DB>",
                CITATION_MMSEQS2, {{"subsetFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDbAndFlat },
                                          {"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"concatdbs",            concatdbs,            &par.concatdbs,            COMMAND_SET,
                "Concatenate two DBs, giving new IDs to entries from 2nd DB",
//                "If exist, the auxillary files: _mapping, source and lookup are also concatenated after IDs update of the 2nd DB",
                "# Download two sequences databases and concat them\n"
                "mmseqs databases PDB pdbDB tmp\n"
                "mmseqs UniProtKB/Swiss-Prot swissprotDB tmp\n"
                "# Works only single threaded since seq. and header DB need the same ordering\n"
                "mmseqs concatdbs pdbDB swissprotDB pdbAndSwissprotDB --threads 1\n"
                "mmseqs concatdbs pdbDB_h swissprotDB_h pdbAndSwissprotDB_h --threads 1\n",
                "Clovis Galiez, Eli Levy Karin & Martin Steinegger (martin.steinegger@snu.ac.kr)",
                "<i:DB> <i:DB> <o:DB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"splitdb",              splitdb,              &par.splitdb,              COMMAND_SET,
                "Split DB into subsets",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:DB> <o:DB>",
                CITATION_MMSEQS2,{{"allDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                         {"allDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"mergedbs",             mergedbs,             &par.mergedbs,             COMMAND_SET,
                "Merge entries from multiple DBs",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:DB> <o:DB> <i:DB1> ... <i:DBn>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::allDb }}},
        {"subtractdbs",          subtractdbs,          &par.subtractdbs,          COMMAND_SET,
                "Remove all entries from first DB occurring in second DB by key",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:resultDBLeft> <i:resultDBRight> <o:resultDB>",
                CITATION_MMSEQS2, {{"resultDBLeft", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"resultDBRight", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},



        {"view",                 view,                 &par.view,                 COMMAND_DB,
                "Print DB entries given in --id-list to stdout",
                "# Print entries with keys 1, 2 and 3 from a sequence DB to stdout\n"
                "mmseqs view sequenecDB --id-list 1,2,3\n",
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:DB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"apply",                apply,                &par.threadsandcompression,
#ifdef __CYGWIN__
                COMMAND_HIDDEN,
#else
                COMMAND_DB,
#endif
                "Execute given program on each DB entry",
                "# Gather all sequences from a cluster DB\n"
                "mmseqs createseqfiledb sequenceDB clusterDB unalignedDB --min-sequences 2\n"
                "# Build MSAs with Clustal-Omega\n"
                "mmseqs apply unalignedDB msaDB -- clustalo -i - -o stdout --threads=1\n\n"
                "# Count lines in each DB entry inefficiently (result2stats is way faster)\n"
                "mmseqs apply DB wcDB -- awk '{ counter++; } END { print counter; }'\n",
                "Milot Mirdita <milot@mirdita.de>",
                "<i:DB> <o:DB> -- program [args...]",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                                           {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"filterdb",             filterdb,             &par.filterDb,             COMMAND_DB,
                "DB filtering by given conditions",
                "# Retain top alignment for each query (alignment DBs are sorted by E-value)\n"
                "mmseqs filterdb alignmentDB topHitAlignmentDB --extract-lines 1\n\n"
                "# Extract alignments with Seq.id. greater than 90%\n"
                "mmseqs filterdb alignmentDB scoreGreater35AlignmentDB --comparison-operator ge --comparison-value 0.9 --filter-column 2\n\n"
                "# Retain all hits matching a regular expression\n"
                "mmseqs filterdb alignmentDB regexFilteredDB --filter-regex '^[1-9].$' --filter-column 2\n\n"
                "# Remove all hits to target keys contained in file db.index\n"
                "mmseqs filterdb --filter-file db.index --positive-filter false\n\n"
                "# Retain all hits matching any boolean expression\n"
                "mmseqs filterdb --filter-expression '$1 * $2 >= 200'\n",
                "Clovis Galiez & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"swapdb",               swapdb,               &par.swapdb,               COMMAND_DB,
                "Transpose DB with integer values in first column",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>, Clovis Galiez & Eli Levy Karin",
                "<i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"prefixid",             prefixid,             &par.prefixid,             COMMAND_DB,
                "For each entry in a DB prepend the entry key to the entry itself",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:DB> <o:DB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"suffixid",             suffixid,             &par.prefixid,             COMMAND_DB,
                "For each entry in a DB append the entry key to the entry itself",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},
        {"renamedbkeys",         renamedbkeys,         &par.renamedbkeys,         COMMAND_DB,
                "Create a new DB with original keys renamed",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:idMapFile|stdin> <i:DB> <o:DB>",
                CITATION_MMSEQS2, {{"idMapFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfileAndStdin },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::allDb }}},

        {"extractorfs",          extractorfs,          &par.extractorfs,          COMMAND_SEQUENCE,
                "Six-frame extraction of open reading frames",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"extractframes",        extractframes,        &par.extractframes,        COMMAND_SEQUENCE,
                "Extract frames from a nucleotide sequence DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> ",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        //TODO remove later?
        {"orftocontig",          orftocontig,          &par.orftocontig,          COMMAND_SEQUENCE,
                "Write ORF locations in alignment format",
                NULL,
                "Eli Levy Karin <eli.levy.karin@gmail.com> ",
                "<i:contigsSequenceDB> <i:extractedOrfsHeadersDB> <o:orfsAlignedToContigDB>",
                CITATION_MMSEQS2, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"reverseseq",          reverseseq,            &par.reverseseq,           COMMAND_SEQUENCE,
                "Reverse (without complement) sequences",
                NULL,
                "Eli Levy Karin <eli.levy.karin@gmail.com> ",
                "<i:sequenceDB> <o:revSequenceDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"translatenucs",        translatenucs,        &par.translatenucs,        COMMAND_SEQUENCE,
                "Translate nucleotides to proteins",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::nuclDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::aaDb }}},
        {"translateaa",          translateaa,          &par.threadsandcompression,COMMAND_SEQUENCE,
                "Translate proteins to lexicographically lowest codons",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::aaDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::nuclDb }}},
        // TODO add SEQUENCE_SPLIT_MODE_SOFT
        {"splitsequence",       splitsequence,         &par.splitsequence,        COMMAND_SEQUENCE,
                "Split sequences by length",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"masksequence",        masksequence,          &par.threadsandcompression,COMMAND_SEQUENCE,
                "Soft mask sequence DB using tantan",
//                "Low. complex regions are masked as lower case characters. The remaining regions are printed as upper case characters.",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"extractalignedregion", extractalignedregion, &par.extractalignedregion, COMMAND_SEQUENCE,
                "Extract aligned sequence region from query",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb}}},



        {"swapresults",          swapresults,          &par.swapresult,           COMMAND_RESULT,
                "Transpose prefilter/alignment DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>, Clovis Galiez & Eli Levy Karin",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::prefAlnResDb },
                                                           {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::prefAlnResDb }}},
        {"result2rbh",           result2rbh,           &par.threadsandcompression,COMMAND_RESULT,
                "Filter a merged result DB to retain only reciprocal best hits",
                NULL,
                "Eli Levy Karin",
                "<i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
        {"result2msa",           result2msa,           &par.result2msa,           COMMAND_RESULT,
                "Compute MSA DB from a result DB",
                NULL,
                "Martin Steinegger (martin.steinegger@snu.ac.kr) & Milot Mirdita <milot@mirdita.de> & Clovis Galiez",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:msaDB>",
                CITATION_MMSEQS2,{{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"msaDB",    DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::msaDb }}},
        {"result2dnamsa",           result2dnamsa,           &par.result2dnamsa,           COMMAND_RESULT,
                "Compute MSA DB with out insertions in the query for DNA sequences",
                NULL,
                "Martin Steinegger (martin.steinegger@snu.ac.kr)",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:msaDB>",
                CITATION_MMSEQS2,{{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                         {"msaDB",    DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::msaDb }}},
        {"result2stats",         result2stats,         &par.result2stats,         COMMAND_RESULT,
                "Compute statistics for each entry in a DB",
                NULL,
                "Clovis Galiez & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:statsDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::allDb },
                                          {"statsDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb }}},
        {"filterresult",         filterresult,         &par.filterresult,         COMMAND_RESULT,
                "Pairwise alignment result filter",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:resultDB>",
                CITATION_MMSEQS2,{{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                         {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
        {"offsetalignment",      offsetalignment,      &par.offsetalignment,      COMMAND_RESULT,
                "Offset alignment by ORF start position",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:queryOrfDB> <i:targetDB> <i:targetOrfDB> <i:alnDB> <o:alnDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"queryOrfDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetOrfDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"alnDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"alnDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"proteinaln2nucl",      proteinaln2nucl,      &par.proteinaln2nucl,      COMMAND_RESULT,
                "Transform protein alignments to nucleotide alignments",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr> ",
                "<i:nuclQueryDB> <i:nuclTargetDB> <i:aaQueryDB> <i:aaTargetDB> <i:alnDB> <o:alnDB>",
                CITATION_MMSEQS2, {{"nuclQueryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"nuclTargetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"aaQueryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"aaTargetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"alnDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"alnDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"result2repseq",       result2repseq,         &par.result2repseq,        COMMAND_RESULT,
                "Get representative sequences from result DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> <i:resultDB> <o:sequenceDb>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"sequenceDb", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"sortresult",           sortresult,           &par.sortresult,           COMMAND_RESULT,
                "Sort a result DB in the same order as the prefilter or align module",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:resultbDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
        {"summarizealis",      summarizealis,      &par.threadsandcompression,     COMMAND_RESULT,
                "Summarize alignment result to one row (uniq. cov., cov., avg. seq. id.)",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:alignmentDB> <o:summerizedDB>",
                CITATION_MMSEQS2|CITATION_UNICLUST, {{"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"summerizedDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb }}},
        {"summarizeresult",      summarizeresult,      &par.summarizeresult,      COMMAND_RESULT,
                "Extract annotations from alignment DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:alignmentDB> <o:alignmentDB>",
                CITATION_MMSEQS2|CITATION_UNICLUST, {{"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},



        {"result2profile",       result2profile,       &par.result2profile,       COMMAND_PROFILE,
                "Compute profile DB from a result DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:profileDB>",
                CITATION_MMSEQS2,{{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"profileDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::profileDb }}},
        {"msa2result",          msa2result,            &par.msa2profile,          COMMAND_PROFILE | COMMAND_EXPERT,
                "Convert a MSA DB to a profile DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:msaDB> <o:seqDB> <o:profileDB>",
                CITATION_SERVER |CITATION_MMSEQS2, {{"msaDB",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::msaDb },
                                                           {"seqDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::aaDb },
                                                           {"profileDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::profileDb }}},
        {"msa2profile",          msa2profile,          &par.msa2profile,          COMMAND_PROFILE | COMMAND_DATABASE_CREATION,
                "Convert a MSA DB to a profile DB",
                "# Convert globally aligned MSAs to profiles\n"
                "# Defines columns as match columns if more than 50% of residues are not gaps\n"
                "# Non-match columns are discarded\n"
                "mmseqs msa2profile msaDB profileDB --match-mode 1 --match-ratio 0.5\n\n"
                "# Assign match-columns through the first sequence\n"
                "# Gaps in query sequence define non-match columns and are discarded\n"
                "mmseqs msa2profile msaDB profileDB --match-mode 0\n",
                "Milot Mirdita <milot@mirdita.de>",
                "<i:msaDB> <o:profileDB>",
                CITATION_SERVER |CITATION_MMSEQS2, {{"msaDB",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::msaDb },
                                                           {"profileDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::profileDb }}},
        {"profile2pssm",         profile2pssm,         &par.profile2pssm,         COMMAND_PROFILE,
                "Convert a profile DB to a tab-separated PSSM file",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:profileDB> <o:pssmFile>",
                CITATION_MMSEQS2, {{"profileDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::profileDb },
                                                           {"pssmFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb }}},
        {"profile2consensus",    profile2consensus,    &par.profile2seq,          COMMAND_PROFILE,
                "Extract consensus sequence DB from a profile DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:profileDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"profileDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::profileDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::aaDb }}},
        {"profile2repseq",       profile2repseq,       &par.profile2seq,          COMMAND_PROFILE,
                "Extract representative sequence DB from a profile DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:profileDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"profileDB",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::profileDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::aaDb }}},
        {"convertprofiledb",     convertprofiledb,     &par.convertprofiledb,     COMMAND_PROFILE,
                "Convert a HH-suite HHM DB to a profile DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:hhsuiteHHMDB> <o:profileDB>",
                CITATION_MMSEQS2,{{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},



        {"enrich",                enrich,              &par.enrichworkflow,       COMMAND_PROFILE_PROFILE,
                "Boost diversity of search result",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                          {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"result2pp",            result2pp,            &par.result2pp,            COMMAND_PROFILE_PROFILE,
                "Merge two profile DBs by shared hits",
                NULL,
                "Clovis Galiez & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <o:profileDB>",
                CITATION_MMSEQS2,{{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::profileDb },
                                         {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::profileDb },
                                         {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                         {"profileDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::profileDb }}},
        {"profile2cs",         profile2cs,             &par.profile2cs,           COMMAND_PROFILE_PROFILE,
                "Convert a profile DB into a column state sequence DB",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:profileDB> <o:csDB>",
                CITATION_MMSEQS2, {{"profileDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::profileDb },
                                         {"csDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::csDb }}},
        {"convertca3m",          convertca3m,          &par.threadsandcompression,COMMAND_PROFILE_PROFILE,
                "Convert a cA3M DB to a result DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:ca3mDB> <o:resultDB>",
                CITATION_MMSEQS2, {{"ca3mDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::ca3mDb },
                                          {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
        {"expandaln",           expandaln,             &par.expandaln,            COMMAND_PROFILE_PROFILE,
                "Expand an alignment result based on another",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <i:resultDB|ca3mDB> <o:alignmentDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"alignmentDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"expand2profile",      expand2profile,        &par.expand2profile,       COMMAND_PROFILE_PROFILE,
                "Expand an alignment result based on another and create a profile",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:queryDB> <i:targetDB> <i:resultDB> <i:resultDB|ca3mDB> <o:profileDB>",
                CITATION_MMSEQS2, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                          {"profileDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::profileDb }}},


        {"diffseqdbs",           diffseqdbs,           &par.diff,                 COMMAND_SPECIAL,
                "Compute diff of two sequence DBs",
//                "It creates 3 filtering files, that can be used in conjunction with \"createsubdb\" tool.\nThe first file contains the keys that has been removed from DBold to DBnew.\nThe second file maps the keys of the kept sequences from DBold to DBnew.\nThe third file contains the keys of the sequences that have been added in DBnew.",
                NULL,
                "Clovis Galiez & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:oldSequenceDB> <i:newSequenceDB> <o:rmSeqKeysFile> <o:keptSeqKeysFile> <o:newSeqKeysFile>",
                CITATION_MMSEQS2, {{"oldSequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"newSequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"rmSeqKeysFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"keptSeqKeysFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"newSeqKeysFile", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"summarizetabs",        summarizetabs,        &par.summarizetabs,        COMMAND_SPECIAL,
                "Extract annotations from HHblits BLAST-tab-formatted results",
                NULL,
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:blastTabDB> <i:lengthFile> <o:summarizedBlastTabDB>",
                CITATION_MMSEQS2|CITATION_UNICLUST,{{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"gff2db",               gff2db,               &par.gff2db,               COMMAND_SPECIAL,
                "Extract regions from a sequence database based on a GFF3 file",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:gff3File1> ... <i:gff3FileN> <i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"gff3File", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                                           {"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"maskbygff",            maskbygff,            &par.gff2db,               COMMAND_SPECIAL,
                "Mask out sequence regions in a sequence DB by features selected from a GFF3 file",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:gff3File> <i:sequenceDB> <o:sequenceDB>",
                CITATION_MMSEQS2, {{"gff3File", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                                           {"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
        {"convertkb",            convertkb,            &par.convertkb,            COMMAND_SPECIAL,
                "Convert UniProtKB data to a DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:uniprotkb.dat[.gz]> ... <i:uniprotkb.dat[.gz]> <o:uniprotkbDB>",
                CITATION_MMSEQS2, {{"DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile },
                                                           {"DB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb }}},
        {"summarizeheaders",     summarizeheaders,     &par.summarizeheaders,     COMMAND_SPECIAL,
                "Summarize FASTA headers of result DB",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:queryDB> <i:targetDB> <i:resultDb> <o:headerDB>",
                CITATION_MMSEQS2|CITATION_UNICLUST, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDb", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"headerDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb}}},
        {"nrtotaxmapping",       nrtotaxmapping,       &par.onlythreads,          COMMAND_SPECIAL,
                "Create taxonomy mapping for NR database",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:accession2taxid1> ... <i:accession2taxidN> <i:seqDB> <o:tsvFile>",
                CITATION_TAXONOMY, {{"ncbiTaxDir", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::VARIADIC, &DbValidator::flatfileAndStdin },
                                    {"nrSeqDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                    {"mappingTSV", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile }}},
        {"extractdomains",       extractdomains,       &par.extractdomains,       COMMAND_SPECIAL,
                "Extract highest scoring alignment regions for each sequence from BLAST-tab file",
                NULL,
                "Milot Mirdita <milot@mirdita.de> & Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:alignmentDb> <i:msaDB> <o:domainDB>",
                CITATION_MMSEQS2|CITATION_UNICLUST, {{"alignmentDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                                                           {"msaDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::msaDb },
                                                           {"domainDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb}}},
        {"countkmer",           countkmer,             &par.countkmer,            COMMAND_SPECIAL,
                "Count k-mers",
                NULL,
                "Martin Steinegger <martin.steinegger@snu.ac.kr>",
                "<i:sequenceDB> ",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},



        {"dbtype",              dbtype,                &par.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_MMSEQS2, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"version",              versionstring,        &par.empty,                COMMAND_HIDDEN,
                "",
                NULL,
                "",
                "",
                CITATION_MMSEQS2, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}},
        {"diskspaceavail",       diskspaceavail,       &par.empty,                COMMAND_HIDDEN,
                "Show available disk space in bytes",
                NULL,
                "",
                "",
                CITATION_MMSEQS2, {{"",DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}}}
};
