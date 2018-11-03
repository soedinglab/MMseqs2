# MMseqs2: ultra fast and sensitive protein search and clustering suite
MMseqs2 (Many-against-Many sequence searching) is a software suite to search and cluster huge protein sequence sets. MMseqs2 is open source GPL-licensed software implemented in C++ for Linux, MacOS, and (as beta version, via cygwin) Windows. The software is designed to run on multiple cores and servers and exhibits very good scalability. MMseqs2 can run 10000 times faster than BLAST. At 100 times its speed it achieves almost the same sensitivity. It can perform profile searches with the same sensitivity as PSI-BLAST at over 400 times its speed.

The MMseqs2 user guide is available in our [GitHub Wiki](https://github.com/soedinglab/mmseqs2/wiki) or as a [PDF file](https://mmseqs.com/latest/userguide.pdf) (Thanks to [pandoc](https://github.com/jgm/pandoc)!)

Please cite: 

[Steinegger M and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi: 10.1038/nbt.3988 (2017)](https://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3988.html).

[Steinegger M and Soeding J. Clustering huge protein sequence sets in linear time. Nature Communications, doi: 10.1038/s41467-018-04964-5 (2018)](https://www.nature.com/articles/s41467-018-04964-5).

![alt tag](https://codeship.com/projects/58db4570-5f19-0134-0f23-2e28d2b4319e/status?branch=master)
![alt tag](https://ci.appveyor.com/api/projects/status/lq8nxeb0j8v38d1a?svg=true)
![alt tag](https://travis-ci.org/soedinglab/MMseqs2.svg?branch=master)
![alt tag](https://zenodo.org/badge/DOI/10.5281/zenodo.840208.svg)

<p align="center"><img src="https://raw.githubusercontent.com/soedinglab/mmseqs2/master/mmseqs2_logo.png" height="256" /></p>

## News
Keep posted about MMseqs2/Linclust updates by following Martin on [Twitter](https://twitter.com/thesteinegger).

07/07/2018 Linclust has just been published at [Nature Communications](https://www.nature.com/articles/s41467-018-04964-5).

17/10/2017 MMseqs2 has just been published at [Nature Biotechnology](https://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3988.html).

19/12/2016 MMseqs2 has a mascot now. "Little Marv" was lovingly crafted by Yuna Kwon. Thank you so much.

## Installation
MMseqs2 can be used by compiling from source, downloading a statically compiled version, using [Homebrew](https://github.com/Homebrew/brew) or [Docker](https://github.com/moby/moby). MMseqs2 requires a 64-bit system (check with `uname -a | grep x86_64`) with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on MacOS).

### Compile from source
Compiling MMseqs2 from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile MMseqs2 `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are needed. Afterwards, the MMseqs2 binary will be located in the `build/bin/` directory.

        git clone https://github.com/soedinglab/MMseqs2.git
        cd MMseqs2
        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
        make
        make install 
        export PATH=$(pwd)/bin/:$PATH

:exclamation: To compile MMseqs2 on MacOS, first install the `gcc` compiler from Homebrew. The default MacOS `clang` compiler does not support OpenMP and MMseqs2 will only be able to use a single thread. Then use the following cmake call:

        CXX="$(brew --prefix)/bin/g++-8" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
                
### Install static Linux version
The following command will download the lastest MMseqs2 Linux version, extract it and set the `PATH` variable.

If your computer supports AVX2 use this (faster than SSE4.1, check by executing `cat /proc/cpuinfo | grep avx2` on Linux and `sysctl -a | grep machdep.cpu.leaf7_features | grep AVX2` on MacOS):

        wget https://mmseqs.com/latest/mmseqs-static_avx2.tar.gz 
        tar xvzf mmseqs-static_avx2.tar.gz
        export PATH=$(pwd)/mmseqs2/bin/:$PATH
        
If your computer supports SSE4.1 use:

        wget https://mmseqs.com/latest/mmseqs-static_sse41.tar.gz 
        tar xvzf mmseqs-static_sse41.tar.gz
        export PATH=$(pwd)/mmseqs2/bin/:$PATH

MMseqs2 comes with a bash command and parameter auto completion, which can be activated by adding the following lines to your $HOME/.bash_profile:

<pre>
        if [ -f /<b>Path to MMseqs2</b>/util/bash-completion.sh ]; then
            source /<b>Path to MMseqs2</b>/util/bash-completion.sh
        fi
</pre>

We also provide static binaries for MacOS and Windows at [mmseqs.com/latest](https://mmseqs.com/latest).

### Install with Homebrew
You can install the latest stable MMseqs2 version for MacOS through [Homebrew](https://github.com/Homebrew/brew) or [Linuxbrew](https://github.com/Linuxbrew/brew) by executing the following:

        brew install mmseqs2

Or the latest development version with:

        brew install https://raw.githubusercontent.com/soedinglab/mmseqs2/master/Formula/mmseqs2.rb --HEAD

This will also automatically install the bash completion (you might have to execute `brew install bash-completion` first).

### Use the Docker image
You can either pull the official docker image by running:

        docker pull soedinglab/mmseqs2

Or build the docker image from the git repository by executing:
        
        git clone https://github.com/soedinglab/MMseqs2.git
        cd MMseqs2
        docker build -t mmseqs2 .

## How to search
You can use the query database "QUERY.fasta" and target database "DB.fasta" in the examples folder to test the search workflow. First, you need to convert the FASTA files into the MMseqs2 database format.

        mmseqs createdb examples/QUERY.fasta queryDB
        mmseqs createdb examples/DB.fasta targetDB
        
If the target database will be used several times, we recommend to precompute an index of `targetDB` as this saves overhead computations. The index should be created on a computer that has the at least the same amount of memory as the computer that performs the search.

Transfer of large database files via NFS quickly becomes time-limiting for MMseqs2. Therefore, ideally the database and database index file should be stored on a fast local drive.

        mmseqs createindex targetDB tmp

MMseqs2 will create, if it does not exist already, a temporary directory `tmp` in which intermediate results are stored. You can also specify a different path, for example on a local drive to reduce load on a shared filesystem or to provide a fast local drive.

:exclamation: In MPI mode all databases and temporary directory need to be accessible by all nodes.

The `mmseqs search` searches the `queryDB` against the `targetDB`. The sensitivity can be adjusted with `-s` parameter and should be adapted based on your use case (see [setting sensitivity -s parameter](https://github.com/soedinglab/mmseqs2/wiki#set-sensitivity--s-parameter)). If you require the exact alignment information in later steps add the option `-a`, without this parameter MMseqs2 will automatically decide if the exact alignment boundaries need to be saved to disk.

Please ensure that, in case of large input databases, the `tmp` directory provides enough free space.
Our user guide provides or information about [disk space requirements](https://github.com/soedinglab/mmseqs2/wiki#prefiltering-module).

To run the search execute:

        mmseqs search queryDB targetDB resultDB tmp

Then convert the result database into a BLAST-tab formatted database (format: qId, tId, seqIdentity, alnLen, mismatchCnt, gapOpenCnt, qStart, qEnd, tStart, tEnd, eVal, bitScore).

        mmseqs convertalis queryDB targetDB resultDB resultDB.m8

Use the option `--format-output "query,target,qaln,taln"` to return query and target accession and the pairwise alignments in tab separated format. You can choose many different [output columns](https://github.com/soedinglab/mmseqs2/wiki#custom-alignment-format-with-convertalis) in the `convertalis` module. Make sure that you used the option `-a` to search (`mmseqs search ... -a`).

        mmseqs convertalis queryDB targetDB resultDB resultDB.pair --format-output "query,target,qaln,taln"

### Other search modes

MMseqs2 provides many additional search modes:
 * Iterative sequences-profile searches (like PSI-BLAST) with the `--num-iterations` parameter
 * [Translated searches](https://github.com/soedinglab/MMseqs2/wiki#translated-sequence-searching) of nucleotides against proteins or proteins against nucleotides
 * [Iterative increasing sensitivity searches](https://github.com/soedinglab/MMseqs2/wiki#how-to-find-the-best-hit-the-fastest-way) to find only the best hits.
 * Fast ungapped alignment searches to find [very similar sequence matches](https://github.com/soedinglab/MMseqs2/wiki#mapping-very-similar-sequences-using-mmseqs-map)
 * Searches against [profile databases such as the PFAM](https://github.com/soedinglab/MMseqs2/wiki#how-to-create-a-target-profile-database-from-pfam)

Many modes can also be combined. You can, for example, do a translated nucleotide against protein profile search.

## How to cluster 
Before clustering, convert your database into the MMseqs2 database format:

        mmseqs createdb examples/DB.fasta DB

Then execute the clustering:

        mmseqs cluster DB clu tmp

Please ensure that in case of large input databases the temporary direcotry provides enough free space.
For disk space requirements, see the user guide.

To generate a FASTA-style formatted output file from the ffindex output file, type:

        mmseqs createseqfiledb DB clu clu_seq 
        mmseqs result2flat DB DB clu_seq clu_seq.fasta
        
To generate a TSV-style formatted output file from the ffindex output file, type:

        mmseqs createtsv DB DB clu clu.tsv
        
To extract the representative sequences from the clustering result call:    
    
        mmseqs result2repseq DB clu DB_clu_rep
        mmseqs result2flat DB DB DB_clu_rep DB_clu_rep.fasta --use-fasta-header

Read more about the format [here](https://github.com/soedinglab/mmseqs2/wiki#clustering-format).

### Memory Requirements
When using MMseqs2 the available memory limits the size of database you will be able to compute in one go.
We recommend at least 128 GB of RAM so you can compute databases up to 30.000.000 entries.
MMseqs2 will automatically subdivide the target database if less memory is available. Runtimes will slightly increase in this case.

You can calculate the memory requirements in bytes for `L` columns and `N` rows using the following formula:
        
        M = (7 × N × L) byte + (8 × a^k) byte

MMseqs2 stores an index table and two auxiliary arrays, which have a total size of `M byte`.

For a database containing `N` sequences with an average length `L`, the memory consumption of the index table is `(7 × N × L) byte` .
Note that the memory consumption grows linearly with the number of the sequences `N` in the database.

The two auxiliary arrays consume `(8 × a^k) byte`, with `a` being the size of the amino acid alphabet (usually 20, the unknown amino acid X is excluded) and the k-mer size `k`.

### How to run MMseqs2 on multiple servers using MPI
MMseqs2 can run on multiple cores and servers using OpenMP and Message Passing Interface (MPI).
MPI assigns database splits to each compute node and each node computes them using multiple cores (OpenMP).
Most of the resource demanding modules of MMseqs2 such as `prefilter` or `align` can take advantage of MPI to speed up the computation.

To parallelize the time-consuming k-mer matching and gapless alignment stages `prefilter` among multiple servers, two different modes are available. In the first, MMseqs2 can split the target sequence set into approximately equal-sized chunks, and each server searches all queries against its chunk. Alternatively, the query sequence set is split into equal-sized chunks and each server searches its query chunk against the entire target set. Splitting the target database is less time-efficient due to the slow, IO-limited merging of results. But it reduces the memory required on each server to `(7 × N L/#chunks) byte + (a^k × 8) byte` and allows users to search through huge databases on servers with moderate memory sizes. If the number of chunks is larger than the number of servers, chunks will be distributed among servers and processed sequentially. By default, MMseqs2 automatically decides which mode to pick based on the available memory (assume that all machines have the same amount of memory). 

Make sure that MMseqs2 was compiled with MPI by using the `-DHAVE_MPI=1` flag (`cmake -DHAVE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..`). Our precompiled static version of MMseqs2 can not use MPI.

To search with multiple server call the `search` workflow with the MPI command exported in the RUNNER environment variable. The databases and temporary folder have to be shared between all nodes (e.g. through NFS):

        RUNNER="mpirun -np 42" mmseqs search queryDB targetDB resultDB tmp

The same requirements apply to clustering or any of the other workflows:

        RUNNER="mpirun -np 42" mmseqs cluster DB clu tmp
