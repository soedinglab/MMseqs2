# MMseqs2.0: ultra fast and sensitive protein search and clustering suite
MMseqs2 (Many-against-Many searching) is a software suite to search and cluster huge protein sequence sets. MMseqs2 is open source GPL-licensed software implemented in C++ for Linux and MacOS. The software is designed to run on multiple cores and servers and exhibits very good scalability. MMseqs2 can run 10000 times faster than BLAST. At 100 times its speed it achieves the same sensitivity. It can also perform profile searches with the same sensitivity as PSI-BLAST but at around 270 times its speed.

Please cite: [Steinegger M and Soeding J. Sensitive protein sequence searching for the analysis of massive data sets. bioRxiv, doi: 10.1101/079681 (2017)](http://biorxiv.org/content/early/2017/03/24/079681).

![alt tag](https://codeship.com/projects/58db4570-5f19-0134-0f23-2e28d2b4319e/status?branch=master)
![alt tag](https://ci.appveyor.com/api/projects/status/lq8nxeb0j8v38d1a?svg=true)
![alt tag](https://travis-ci.org/soedinglab/MMseqs2.svg?branch=master)
![alt tag](https://zenodo.org/badge/DOI/10.5281/zenodo.840208.svg)

<p align="center"><img src="https://raw.githubusercontent.com/soedinglab/mmseqs2/master/mmseqs2_logo.png" height="256" /></p>


## News
Keep posted about MMseqs2/Linclust updates by following Martin on [twitter](https://twitter.com/thesteinegger).

04/08/2017 We had a serious bug in the linclust/cluster workflow when performing target-cov based clustering since 20/07/2017. We will also add the clustering workflows to our regression tests. The bug should be fixed now. Please update your MMseqs2.

05/25/2017 We updated the Linclust manuscript. Linclust is now 3x faster and we added an metagenomic protein assembly application. Happy towel day. A pre-print can be downloaded here: [Steinegger M and Soeding J. Linclust: clustering billions of protein sequences per day on a single server (2017)](http://biorxiv.org/content/early/2017/05/25/104034). 

30/01/2017 We added a new clustering workflow called "Linclust". Linclust can cluster sequences in linear time down to 50% sequence identity. The Metaclust95 and Metaclust50 database can be download at [metaclust.mmseqs.com](https://metaclust.mmseqs.com/). A pre-print can be downloaded here: [Steinegger M and Soeding J. Linclust: clustering protein sequences in linear time (2017)](http://www.biorxiv.org/content/early/2017/01/29/104034.article-metrics). 

19/12/2016 MMseqs2 has a mascot now. It is the "little Marv" and was lovingly crafted by Yuna Kwon. Thank you so much.

07/12/2016 We added a new parameter called --max-accept. This parameter limits the amount of alignments that get accepted. Please do not use --max-seqs to limit your result size since it decreases the sensitivity of MMseqs2.

## Documentation
The MMseqs2 user guide is available in our [Github Wiki](https://github.com/soedinglab/mmseqs2/wiki). A PDF version of the user guide is also [available](https://mmseqs.com/latest/userguide.pdf) (Thanks to [pandoc](https://github.com/jgm/pandoc)!).
You can read on for a quick start guide with installation instructions and examples for searching and clustering.

## Installation
MMseqs can be installed by compiling the binary from source, download a statically compiled version, using [Homebrew](https://github.com/Homebrew/brew) or [Docker](https://github.com/moby/moby). MMseqs2 requires a 64-bit system (check with `uname -a | grep x86_64`) with at least the SSE4.1 intruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux and `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on MacOS).

### Compile from source
Compiling MMseqs2 from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile MMseqs2 `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are needed. Afterwards, the MMseqs2 binary will be located in `build/bin/`.

        git clone https://github.com/soedinglab/MMseqs2.git
        cd MMseqs2
        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
        make
        make install 
        export PATH=$(pwd)/bin/:$PATH
        
:exclamation: Please install and use `gcc` from Homebrew, if you want to compile MMseqs2 on MacOS. The default MacOS `clang` compiler does not support OpenMP and MMseqs2 will not be able to run multithreaded. Use the following cmake call:

        CXX="$(brew --prefix)/bin/g++-6" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
                
### Install static Linux version
The following command will download the last MMseqs version, extract it and set the `PATH` variable. This version runs only on linux. If you want to run it on Mac please compile it or use brew.

If your computer supports AVX2 use this (faster than SSE4.1, check by executing `cat /proc/cpuinfo | grep avx2` on Linux and `sysctl -a | grep machdep.cpu.leaf7_features | grep AVX2` on MacOS):

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

### Install with Homebrew
You can install MMseqs2 for Mac OS through [Homebrew](https://github.com/Homebrew/brew) by executing the following:

        brew install https://raw.githubusercontent.com/soedinglab/mmseqs2/master/Formula/mmseqs.rb --HEAD

This will also automatically install the bash completion (you might have to execute `brew install bash-completion` first). This will also work for [Linuxbrew](https://github.com/Linuxbrew/brew).

### Use the Docker image
You can either pull the official docker image by running:

        docker pull soedinglab/mmseqs2

Or build the docker image from the git repository by executing:
        
        git clone https://github.com/soedinglab/MMseqs2.git
        cd MMseqs2
        docker build -t mmseqs2 .

## How to search
You can use the query database "queryDB.fasta" and target database "targetDB.fasta" in the examples folder to test the search workflow. First, you need to convert the fasta files into mmseqs database format. 

        mmseqs createdb examples/QUERY.fasta queryDB
        mmseqs createdb examples/DB.fasta targetDB
        
If the target database is to be used several times, it is recommended to precompute an index of the targetDB as this saves overhead computations. The index should be created on a computer that has the same amount of memory as the computer that performs the search. 

Transfer of large database files via NFS quickly becomes time-limiting for MMseqs2. Therefore, ideally the database and database index file should be stored on a fast local drive.

        mmseqs createindex targetDB
        
You need to create a temporary directory in which MMseqs2 will store intermediate results.

        mkdir <fast_local_drive>/tmp

It is recommend to create this tmp on a local drive to reduce load on the NFS.

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
        
        M = (7 × N × L) byte + (8 × a^k) byte

MMseqs stores an index table and two auxiliary arrays, which have a total size of `M byte`.

For a database containing N sequences with an average length L, the memory consumption of the index table is `(7 × N × L) byte` .
Note that the memory consumption grows linearly with the number of the sequences N in the database.

The two auxiliary arrays consume `(8 × a^k) byte`, with `a` being the size of the amino acid alphabet (usually 21 including the unknown amino acid X) and the k-mer size `k`.

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

