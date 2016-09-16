# MMseqs2.0: ultra fast and sensitive search and clustering suite
MMseqs2 is a software suite for very fast protein sequence searches and clustering of huge protein sequence data sets. MMseqs2 reaches the same sensitivity as BLAST magnitude faster and which can also perform profile searches like PSI-BLAST but also ~500x faster.

MMseqs2 has not just improved in sensitivity and speed it also improved usability through several helper tools, to convert formats,  do six frame translations, GFF tools, ...

## Installation
The manuscript for MMseqs2 is still in preparations therefore we can just offer static compiled binaries. We hope that the full source code will be available within the next weeks.

### Static version
The following command will download the last MMseqs version, extract it and set the environment variables.

        wget https://mmseqs.com/latest/mmseqs.tar.gz
        tar xvfz mmseqs.tar.gz

MMseqs comes with a bash command and parameter auto completion
by pressing tab. The bash completion for subcommands and parameters can be installed by adding the following lines to your $HOME/.bash_profile:

        if [ -f /path/to/mmseqs/util/bash-completion.sh ]; then
         .  /path/to/mmseqs/util/bash-completion.sh
        fi

#### Mac OS
You can install MMseqs2 for Mac OS through [Homebrew](https://github.com/Homebrew/brew) by executing the following:

        brew install https://raw.githubusercontent.com/soedinglab/mmseqs2/master/Formula/mmseqs.rb --HEAD

This will also automatically install the bash completion (you might have to do `brew install bash-completion` first).

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


## Requirements

To compile from source, you will need:

  * a recent C and C++ compiler (Minimum requirement is GCC 4.6. GCC 4.8 or later is recommended).

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
MMseqs contains one binary `mmseqs`. It contains several commands to execute complete workflows that combines the MMseqs core modules. 
The other three commands execute the single modules which are used by the workflows and are available for advanced users.

### Workflows
* `mmseqs search` Compares all sequences in the query database with all sequences in the target database.
* `mmseqs cluster` Clusters the sequences in the input database by sequence similarity.

### Single modules
* `mmseqs prefilter` Computes k-mer similarity scores between all sequences in the query database and all sequences in the target database.
* `mmseqs align` Computes Smith-Waterman alignment scores between all sequences in the query database and the sequences of the target database whose prefiltering scores
* `mmseqs cluster` Computes a similarity clustering of a sequence database based on Smith-Waterman alignment scores of the sequence pairs computed by `mmseqs align`.

