# MMseqs2.0: ultra fast and sensitive search and clustering suite
MMseqs2 (Many-against-Many searching) is a software suite to search and cluster huge sequence sets. MMseqs2 is open source GPL-licensed software implemented in C++ for Linux and Mac OS. The software is designed to run on multiple cores and servers and exhibits very good scalability. MMseqs2 reaches the same sensitivity as BLAST magnitude faster and which can also perform profile searches like PSI-BLAST but also ~270x faster.

Please cite M. Steinegger, J. Soeding, "Sensitive protein sequence searching for the analysis of massive data sets.", bioRxiv 079681 (2016).

## Installation
MMseqs can be installed by compiling the binary, download a statically compiled version or using [Homebrew](https://github.com/Homebrew/brew). MMseqs2 requires a 64-bit system (check with `uname -a | grep x86_64`) with at least the SSE4.1 intruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux and `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on MacOS).

### Compile
Compiling MMseqs2 from source has the advantage that it will be optimized to the specific system, which might improve its performance. To compile MMseqs2 `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are needed. Afterwards, the MMseqs2 binary will be located in in `build/bin/`.

        git clone https://github.com/soedinglab/MMseqs2.git
        cd MMseqs2
        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
        make
        make install 
        export PATH=$(pwd)/bin/:$PATH
                
### Static Linux version
The following command will download the last MMseqs version, extract it and set the environment variables. This version runs just on linux. If you want to run it on Mac please compile it or use brew.

If your computer supports AVX2 use this (faster than SSE4.1):

        wget https://mmseqs.com/latest/mmseqs-static_avx2.tar.gz 
        tar xvzf mmseqs-static_avx2.tar.gz
        export PATH=$(pwd)/mmseqs/bin/:$PATH
        
If your computer supports SSE4.1 use:

        wget https://mmseqs.com/latest/mmseqs-static_sse41.tar.gz 
        tar xvzf mmseqs-static_sse41.tar.gz
        export PATH=$(pwd)/mmseqs/bin/:$PATH

MMseqs comes with a bash command and parameter auto completion
by pressing tab. The bash completion for subcommands and parameters can be installed by adding the following lines to your $HOME/.bash_profile:

<pre>
        if [ -f /<b>Path to MMseqs2</b>/util/bash-completion.sh ]; then
         .  /<b>Path to MMseqs2</b>/util/bash-completion.sh
        fi
</pre>
#### [Homebrew](https://github.com/Homebrew/brew) 
You can install MMseqs2 for Mac OS through [Homebrew](https://github.com/Homebrew/brew) by executing the following:

        brew install https://raw.githubusercontent.com/soedinglab/mmseqs2/master/Formula/mmseqs.rb --HEAD

This will also automatically install the bash completion (you might have to do `brew install bash-completion` first).
The formula will also work for [Linuxbrew](https://github.com/Linuxbrew/brew).

## How to search
You can use the query database "queryDB.fasta" and target database "targetDB.fasta" in the examples folder to test the search workflow.
As a frist step you need to convert the fasta files to into mmseqs database format. The database should be stored on a local drive.

        mmseqs createdb examples/QUERY.fasta queryDB
        mmseqs createdb examples/DB.fasta targetDB
        
It is recommended to precomputes an index of the the targetDB for fast read in if the target databse is reused for again. The index should be created on a computer that has the same amount of memory as the computer that performs the search. 

        mmseqs createindex targetDB
        
MMseqs can produce a high IO on the file system. It is recommend to create this tmp on a local drive.        

        mkdir tmp

The `mmseqs search` searches the `queryDB` against the `targetDB`. The sensitivity can be adjusted with `-s` and should be adapted based on your use case. If you want to use alignment backtraces in later steps add the option `-a`.  An iterative profile search (like PSI-BLAST) can be trigged with `--num-iterations`. 

Please ensure that in case of large input databases tmp provides enough free space.
For the disc space requirements, see the user guide.
To run the search type:

        mmseqs search queryDB targetDB resultDB tmp

Then convert the result database into a BLAST-tab formatted database (format: qId, tId, seqIdentity, alnLen, mismatchCnt, gapOpenCnt, qStart, qEnd, tStart, tEnd, eVal, bitScore).

        mmseqs convertalis queryDB targetDB resultDB resultDB.m8

Use the option `--format-mode 1` to convert the results to pairwise alignments. Make sure that you searched with the option `-a` (`mmseqs search ... -a`).

        mmseqs convertalis queryDB targetDB resultDB resultDB.pair --format-mode 1

## How to cluster 
Before clustering, convert your database into the mmseqs database format:

        mmseqs createdb examples/DB.fasta DB

Create a temporary folder for the run. 

        mkdir tmp

Please ensure that in case of large input databases the temporary folder tmp provides enough free space.
For the disc space requirements, see the user guide. 

        mmseqs cluster DB clu tmp

To generate a FASTA-style formatted output file from the ffindex output file, type:

        mmseqs createseqfiledb DB clu clu_seq 
        mmseqs result2flat DB DB clu_seq clu_seq.fasta
        
To generate a TSV-style formatted output file from the ffindex output file, type:

        mmseqs createtsv DB DB clu clu.tsv

### Memory Requirements
When using MMseqs the available memory limits the size of database you will be able to compute. 
We recommend at least 128 GB of RAM so you can compute databases up to 30.000.000 entries:

You can calculate the memory requirements in bytes for L columns and N rows using the following formula:
        
        M = (7 B × N × L + 8 B × a^k) byte

MMseqs stores an index table and two auxiliary arrays, which have a total size of M byte.

For a database containing N sequences with an average length L, the memory consumption of the index table is `(7B × N × L)` byte.
Note that the memory consumption grows linearly with the number of the sequences N in the database.

The two auxiliary arrays consume `(8 B × a^k)` bytes, with a being the size of the amino acid alphabet (usually 21 including the unknown amino acid X) and the  k-mer size k.

### How to run MMseqs2 on multipe servers using MPI
MMseqs2 can run on multiple cores and servers using OpenMP (OMP) and message passing interface (MPI).
MPI assigns database splits to each servers and each server computes them using multiple cores (OMP). 
Currently `prefilter`, `align`, `result2profile`, `swapresults` can take advantage of MPI.

To parallelize the time-consuming k-mer matching and gapless alignment stages `prefilter` among multiple servers, two different modes are available. In the first, MMseqs2 can split the target sequence set into approximately equal-sized chunks, and each server searches all queries against its chunk. Alternatively, the query sequence set is split into equal-sized chunks and each server searches its query chunk against the entire target set. Splitting the target database is less time-efficient due to the slow, IO-limited merging of results. But it reduces the memory required on each server to `7 × N L/#chunks + a^k × 8 B` and allows users to search through huge databases on servers with moderate memory sizes. If the number of chunks is larger than the number of servers, chunks will be distributed among servers and processed sequentially. By default, MMseqs2 automatically decides which mode to pick based on the available memory (assume that all machines have the same amount of memory). 

Make sure that MMseqs2 was compiled with MPI by using the `-DHAVE_MPI=1` flag (`cmake -DHAVE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..`). Our precomplied static version of MMseqs2 can not use MPI.

To search with multiple server just call the search and add the RUNNER variable. The TMP folder has to be shared between all nodes (e.g. NFS)

        RUNNER="mpirun -np 42" mmseqs search queryDB targetDB resultDB tmp

For clustering just call the clustering. The TMP folder has to be shared between all nodes (e.g. NFS)

        RUNNER="mpirun -np 42" mmseqs cluster DB clu tmp

## Overview of MMseqs
MMseqs2 is a stand-alone binary `mmseqs`, which contains commands to execute complete workflows, tools or utilities. 
This modular architecture, can be used chain tools together to create workflows for analysing huge sequence sets. Three plug-and-play bash-scripted workflows for sequence searching `mmseqs search`, sequence clustering `mmseqs cluster`, and updating clusterings`clusterupdate` facilitate the usage for standard tasks. Example bash scripted workflows can be found in the `data` folder.

### Main tools
* `mmseqs createdb` converts a protein sequence set in a FASTA formatted file to MMseqs’ sequence DB format. This format is needed as input to mmseqs search and many other tools.
* `mmseqs search` Search with query sequence or profile DB (iteratively) through target sequence DB", Searches with the sequences or profiles query DB through the target sequence DB by running the prefilter tool and the align tool for Smith-Waterman alignment. For each query a results file with sequence matches is written as entry into a database of search results (“alignmentDB”). In iterative profile search mode, the detected sequences satisfying user-specified criteria are aligned to the query MSA, and the resulting query profile is used for the next search iteration. Iterative profile searches are usually much more sensitive than (and at least as sensitive as) searches with single query sequences.
* `mmseqs cluster`  Clusters sequences by similarity. It compares all sequences in the sequence DB with each other using mmseqs search, filters alignments according to user-specified criteria (max. E-value, min. coverage,...), and runs mmseqs clust to group similar sequences together into clusters.
* `mmseqs createindex` Precomputes an index table for the sequence DB. Handing over the precomputed index table as input to mmseqs search or mmseqs prefilter eliminates the computational overhead of building the index table on the fly.

### Core modules
* `prefilter` Search with query sequence / profile DB through target DB (k-mer matching + ungapped alignment)
* `align` Compute Smith-Waterman alignments for previous results (e.g. prefilter DB, cluster DB)
* `clust` Cluster sequence DB from alignment DB (e.g. created by searching DB against itself)
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
