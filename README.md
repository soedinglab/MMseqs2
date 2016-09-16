# MMseqs2.0: ultra fast and sensitive search and clustering suite
MMseqs2 (Many-against-Many searching) is a software suite to search and cluster huge sequence sets. MMseqs2 is open source GPL-licensed software implemented in C++ for Linux and Mac OS. The software is designed run on multiple cores and servers and exhibits very good scalability. MMseqs2 reaches the same sensitivity as BLAST magnitude faster and which can also perform profile searches like PSI-BLAST but also ~270x faster.

MMseqs2 has not just improved in sensitivity and speed it also improved usability through several helper tools. The utilities comprise tools for format conversion, multiple sequence alignment, sequence profile calculation, 6-frame translation for ORF extraction, set operations on sequence sets, regex-based filters, and statistics tools to analyse results.

## Installation
### Compile
Compiling MMseqs2 from source has the advantage that it will be optimized to the specific system, which might improve its performance. To compile mmseqs `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are needed. Afterwards, the MMseqs2 binary will be located in in `build/bin/`.

        git clone https://github.com/soedinglab/MMseqs2.git
        cd mmseqs
        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
        make
        make install 
        
        
### Static version
The following command will download the last MMseqs version, extract it and set the environment variables.

        wget https://mmseqs.com/latest/mmseqs.tar.gz
        tar xvzf mmseqs.tar.gz

MMseqs comes with a bash command and parameter auto completion
by pressing tab. The bash completion for subcommands and parameters can be installed by adding the following lines to your $HOME/.bash_profile:

        if [ -f /path/to/mmseqs/util/bash-completion.sh ]; then
         .  /path/to/mmseqs/util/bash-completion.sh
        fi

#### [Homebrew](https://github.com/Homebrew/brew) 
You can install MMseqs2 for Mac OS through [Homebrew](https://github.com/Homebrew/brew) by executing the following:

        brew install https://raw.githubusercontent.com/soedinglab/mmseqs2/master/Formula/mmseqs.rb --HEAD

This will also automatically install the bash completion (you might have to do `brew install bash-completion` first).
The formula will also work for [Linuxbrew](https://github.com/Linuxbrew/brew).

### How to search
You can use the query database queryDB.fasta and target database targetDB.fasta to test the search workflow.
Before clustering, you need to convert your database containing query sequences (queryDB.fasta) and your target database (targetDB.fasta) into mmseqs database format:

        mmseqs createdb example/queryDB.fasta queryDB
        mmseqs createdb example/targetDB.fasta targetDB
        
It generates ffindex database files, e. g. queryDB and ffindex index file queryDB.index
from queryDB.fasta. Then, generate a directory for tmp files:
For the next step computes an index file of the the targetDB for fast read in. It is recommend to compute the index if the targetDB is reused for search several times.

        mmseqs createindex targetDB
        
MMseqs can produce a high IO on the file system. It is recommend to create this tmp on a local drive.        

        mkdir tmp

Please ensure that in case of large input databases tmp provides enough free space.
For the disc space requirements, see the user guide.
To run the search type:

        mmseqs search queryDB targetDB resultDB tmp

Then convert the result ffindex database into a FASTA formatted database: 

        mmseqs convertalis queryDB targetDB resultDB resultDB.m8


### How to cluster 
Before clustering, convert your FASTA database into ffindex format:

        mmseqs createdb $MMDIR/data/DB.fasta DB

Please ensure that in case of large input databases the temporary folder tmp  provides enough free space.
For the disc space requirements, see the user guide. 

        mkdir tmp
        mmseqs cluster DB clu tmp

To generate a FASTA-style formatted output file from the ffindex output file, type:

        mmseqs createseqfiledb DB DB_clu DB_clu_seq 
        mmseqs result2flat DB DB DB_clu_seq DB_clu_seq.fasta
        
To generate a TSV-style formatted output file from the ffindex output file, type:

        mmseqs createtsv DB DB clu clu.tsv


### Memory Requirements
When using MMseqs the available memory limits the size of database you will be able to compute. 
We recommend at least 128 GB of RAM so you can compute databases up to 30.000.000 entries:

You can calculate the memory requirements in bytes for L columns and N rows using the following formula:
        
        M = (7*N*L + 8*a^k) byte

MMseqs stores an index table and two auxiliary arrays, which have a total size of M byte.

For a database containing N sequences with an average length L, the memory consumption of the index table is `(7*N*L) byte`.
Note that the memory consumption grows linearly with the number of the sequences N in the database.

The two auxiliary arrays consume `(8*a^k) byte`, with a being the size of the amino acid alphabet (usually 21 including the unknown amino acid X) and the  k-mer size k.

## Overview of MMseqs
MMseqs is stand-alone binary `mmseqs`, which contains several commands to execute complete workflows, tools or utilities. 
MMseqs modular architecture, can be used chain tools together to create workflows for analysing huge sequence sets. Three plug-and-play bash-scripted workflows for sequence searching `mmseqs search`, sequence clustering `mmseqs cluster`, and updating clusterings`clusterupdate` facilitate the usage for standard tasks. Example bash scripted workflows can be found in the `util` folder.

### Main tools
* `mmseqs createdb` converts a protein sequence set in a FASTA formatted file to MMseqs’ sequence DB format. This format is needed as input to mmseqs search and many other tools.
* `mmseqs search` Search with query sequence or profile DB (iteratively) through target sequence DB", Searches with the sequences or profiles query DB through the target sequence DB by running the prefilter tool and the align tool for Smith-Waterman alignment. For each query a results file with sequence matches is written as entry into a database of search results (“alignmentDB”). In iterative profile search mode, the detected sequences satisfying user-specified criteria are aligned to the query MSA, and the resulting query profile is used for the next search iteration. Iterative profile searches are usually much more sensitive than (and at least as sensitive as) searches with single query sequences.
* `mmseqs cluster`  Clusters sequences by similarity. It compares all sequences in the sequence DB with each other using mmseqs search, filters alignments according to user-specified criteria (max. E-value, min. coverage,...), and runs mmseqs clust to group similar sequences together into clusters.
* `mmseqs createindex` Precomputes an index table for the sequence DB. Handing over the precomputed index table as input to mmseqs search or mmseqs prefilter eliminates the computational overhead of building the index table on the fly.


### Core modules
* `mmseqs prefilter` Search with query sequence / profile DB through target DB (k-mer matching + ungapped alignment)
* `mmseqs align` Compute Smith-Waterman alignments for previous results (e.g. prefilter DB, cluster DB)
* `mmseqs clust` Cluster sequence DB from alignment DB (e.g. created by searching DB against itself)
* `clustlinear`       	Cluster sequences of >70% sequence identity *in linear time*
* `clusthash`         	Cluster sequences of same length and >90% sequence identity *in linear time*

### Utility tools for format conversions
* `createtsv`         	Create tab-separated flat file from prefilter DB, alignment DB, or cluster DB
* `convertalis`       	Convert alignment DB to BLAST-tab format, SAM flat file, or to raw pairwise alignments
* `convertprofiledb`  	Convert ffindex DB of HMM/HMMER3/PSSM files to MMseqs profile DB
* `convert2fasta`     	Convert sequence DB to FASTA format
* `result2flat`       	Create a FASTA-like flat file from prefilter DB, alignment DB, or cluster DB

### Utility tools for clustering
* `clusterupdate`     	Update clustering of old sequence DB to clustering of new sequence DB
* `createseqfiledb`   	Create DB of unaligned FASTA files (1 per cluster) from sequence DB and cluster DB
* `mergeclusters`     	Merge multiple cluster DBs into single cluster DB

### Utility tools to manipulate DBs
* `extractorfs`       	Extract open reading frames from all six frames from nucleotide sequence DB
* `translatenucs`     	Translate nucleotide sequence DB into protein sequence DB
* `swapresults`       	Reformat prefilter/alignment/cluster DB as if target DB had been searched through query DB
* `mergedbs`          	Merge multiple DBs into a single DB, based on IDs (names) of entries
* `splitdb`           	Split a mmseqs DB into multiple DBs
* `subtractdbs`       	Generate a DB with entries of first DB not occurring in second DB
* `filterdb`          	Filter a DB by conditioning (regex, numerical, ...) on one of its whitespace-separated columns
* `createsubdb`       	Create a subset of a DB from a file of IDs of entries
* `result2profile`    	Compute profile and consensus DB from a prefilter, alignment or cluster DB
* `result2msa`        	Generate MSAs for queries by locally aligning their matched targets in prefilter/alignment/cluster DB
* `result2stats`      	Compute statistics for each entry in a sequence, prefilter, alignment or cluster DB

### Special-purpose utilities
* `diffseqdbs`        	Find IDs of sequences kept, added and removed between two versions of sequence DB
* `concatdbs`         	Concatenate two DBs, giving new IDs to entries from second input DB
* `summarizetabs`     	Extract annotations from HHblits BAST-tab-formatted results
* `gff2db`            	Turn a gff3 (generic feature format) file into a gff3 DB
* `maskbygff`         	X out sequence regions in a sequence DB by features in a gff3 file
* `prefixid`          	For each entry in a DB prepend the entry ID to the entry itself
* `convertkb`         	Convert UniProt knowledge flat file into knowledge DB for the selected column types
* `summarizeheaders`  	Return a new summarized header DB from the UniProt headers of a cluster DB
* `extractalignedregion`	Extract aligned sequence region
* `extractdomains`    	Extract highest scoring alignment region for each sequence from BLAST-tab file
