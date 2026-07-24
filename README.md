# MMseqs2: ultra fast and sensitive sequence search and clustering suite
MMseqs2 (Many-against-Many sequence searching) is a software suite to search and cluster huge protein and nucleotide sequence sets. MMseqs2 is free and open source software implemented in C++ for Linux, MacOS, and (as beta version, via cygwin) Windows. The software is designed to run on multiple cores and servers and exhibits very good scalability. MMseqs2 can run 10000 times faster than BLAST. At 100 times its speed it achieves almost the same sensitivity. It can perform profile searches with the same sensitivity as PSI-BLAST at over 400 times its speed.

##  Publications

[Steinegger M and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi: 10.1038/nbt.3988 (2017)](https://www.nature.com/articles/nbt.3988).

[Steinegger M and Soeding J. Clustering huge protein sequence sets in linear time. Nature Communications, doi: 10.1038/s41467-018-04964-5 (2018)](https://www.nature.com/articles/s41467-018-04964-5).

[Mirdita M, Steinegger M and Soeding J. MMseqs2 desktop and local web server app for fast, interactive sequence searches. Bioinformatics, doi: 10.1093/bioinformatics/bty1057 (2019)](https://academic.oup.com/bioinformatics/article/35/16/2856/5280135).

[Mirdita M, Steinegger M, Breitwieser F, Soding J, Levy Karin E: Fast and sensitive taxonomic assignment to metagenomic contigs. Bioinformatics, doi: 10.1093/bioinformatics/btab184 (2021)](https://doi.org/10.1093/bioinformatics/btab184).

[Kallenborn F, Chacon A, Hundt C, Sirelkhatim H, Didi K, Cha S, Dallago C, Mirdita M, Schmidt B, Steinegger M: GPU-accelerated homology search with MMseqs2. Nature Methods, doi: 10.1038/s41592-025-02819-8 (2025)](https://www.nature.com/articles/s41592-025-02819-8)

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/mmseqs2.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/mmseqs2) [![Github All Releases](https://img.shields.io/github/downloads/soedinglab/mmseqs2/total.svg)](https://github.com/soedinglab/mmseqs2/releases/latest) [![Biocontainer Pulls](https://img.shields.io/endpoint?url=https%3A%2F%2Fmmseqs.com%2Fbiocontainer.php%3Fcontainer%3Dmmseqs2)](https://biocontainers.pro/#/tools/mmseqs2) [![Build Status](https://dev.azure.com/themartinsteinegger/mmseqs2/_apis/build/status/soedinglab.MMseqs2?branchName=master)](https://dev.azure.com/themartinsteinegger/mmseqs2/_build/latest?definitionId=2&branchName=master)

<p align="center"><img src="https://raw.githubusercontent.com/soedinglab/mmseqs2/master/.github/mmseqs2_logo.png" height="256" /></p>


## Documentation
The MMseqs2 user guide is available in our [GitHub Wiki](https://github.com/soedinglab/mmseqs2/wiki) or as a [PDF file](https://mmseqs.com/latest/userguide.pdf) (Thanks to [pandoc](https://github.com/jgm/pandoc)!). The wiki also contains [tutorials](https://github.com/soedinglab/MMseqs2/wiki/Tutorials) to learn how to use MMseqs2 with real data. For questions please open an issue on [GitHub](https://github.com/soedinglab/MMseqs2/issues).
Keep posted about MMseqs2 updates by following [Martin](https://bsky.app/profile/martinsteinegger.bsky.social) and [Milot](https://bsky.app/profile/milot.bsky.social).

## Installation
MMseqs2 can be used by [compiling from source](https://github.com/soedinglab/MMseqs2/wiki#installation), downloading a statically compiled binary at [mmseqs.com/latest](https://mmseqs.com/latest), using [Homebrew](https://github.com/Homebrew/brew), [conda](https://github.com/conda-forge/miniforge) or [Docker](https://github.com/moby/moby).

When compiling from source, pass `-DMMSEQS_INT64_IDS=1` to CMake to build MMseqs2 with 64-bit internal sequence/database IDs instead of the default 32-bit IDs.

    # install by brew
    brew install mmseqs2
    # install via conda
    conda install -c conda-forge -c bioconda mmseqs2
    # install docker
    docker pull ghcr.io/soedinglab/mmseqs2
    # MMseqs2-GPU mostly-static AVX2 build requiring glibc >= 2.29 and nvidia driver >=525.60.13 (see below)
    wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz; tar xvfz mmseqs-linux-gpu.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
    # static build with AVX2 (fastest)
    wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
    # static AVX2 build with support for >2^32 sequences. NEED 2x MORE MEMORY EVEN FOR SMALLER DBS!
    wget https://mmseqs.com/latest/mmseqs-linux-avx2-64bit.tar.gz; tar xvfz mmseqs-linux-avx2-64bit.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
    # static build with SSE4.1
    wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz; tar xvfz mmseqs-linux-sse41.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
    # static build with SSE2 (slowest, for very old systems)
    wget https://mmseqs.com/latest/mmseqs-linux-sse2.tar.gz; tar xvfz mmseqs-linux-sse2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
    # macOS (universal, works on Apple Silicon and Intel Macs)
    wget https://mmseqs.com/latest/mmseqs-osx-universal.tar.gz; tar xvfz mmseqs-osx-universal.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

MMseqs2 requires an AMD or Intel 64-bit system (check with `uname -a | grep x86_64`). We recommend using a system with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on MacOS). The AVX2 version is faster than SSE4.1, check if AVX2 is supported by executing `cat /proc/cpuinfo | grep avx2` on Linux and `sysctl -a | grep machdep.cpu.leaf7_features | grep AVX2` on MacOS. A SSE2 version is also available for very old systems. MMseqs2 also works on ARM64 systems and on PPC64LE systems with POWER8 ISA or newer. MMseqs2-GPU requires the support of CUDA-enabled GPUs of the Turing generation or newer.

> [!NOTE]
> We recently added support for GPU-accelerated protein sequence and profile searches. This requires an NVIDIA GPU of the Ampere generation or newer for full speed, however, also works at reduced speed for Turing-generation GPUs. The bioconda- and precompiled binaries will not work on older GPU generations (e.g. Volta or Pascal).
> Check the [wiki](https://github.com/soedinglab/MMseqs2/wiki#compile-from-source-for-linux-with-gpu-support) for instructions on how to get started.

MMseqs2 comes with a bash command and parameter auto completion, which can be activated by adding the following to your $HOME/.bash_profile:

<pre>
if [ -f /<b>Path to MMseqs2</b>/util/bash-completion.sh ]; then
    source /<b>Path to MMseqs2</b>/util/bash-completion.sh
fi
</pre>
         
## Getting started
We provide `easy` workflows to cluster, search and assign taxonomy. These `easy` workflows are a shorthand to deal directly with FASTA/FASTQ files as input and output. MMseqs2 provides many modules to transform, filter, execute external programs and search. However, these modules use the MMseqs2 database formats, instead of the FASTA/FASTQ format. For maximum flexibility, we recommend using MMseqs2 workflows and modules directly. Please read more about this in the [documentation](https://github.com/soedinglab/mmseqs2/wiki).

### Cluster

For clustering, MMseqs2 `easy-cluster` and `easy-linclust` are available.

`easy-cluster` by default clusters the entries of a FASTA/FASTQ file using a cascaded clustering algorithm.
        
    mmseqs easy-cluster examples/DB.fasta clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1
        
`easy-linclust` clusters the entries of a FASTA/FASTQ file. The runtime scales linearly with input size. This mode is recommended for huge datasets.
                
    mmseqs easy-linclust examples/DB.fasta clusterRes tmp
                
Read more about the [clustering format](https://github.com/soedinglab/mmseqs2/wiki#clustering-format) in our user guide.
                
Please adjust the [clustering criteria](https://github.com/soedinglab/MMseqs2/wiki#clustering-criteria) and check if temporary directory provides enough free space. For disk space requirements, see the user guide.

### Search
         
The `easy-search` workflow searches directly with a FASTA/FASTQ files against either another FASTA/FASTQ file or an already existing MMseqs2 database.
        
    mmseqs easy-search examples/QUERY.fasta examples/DB.fasta alnRes.m8 tmp

 
It is also possible to pre-compute the index for the target database. This reduces overhead when searching repeatedly against the same database.

    mmseqs createdb examples/DB.fasta targetDB
    mmseqs createindex targetDB tmp
    mmseqs easy-search examples/QUERY.fasta targetDB alnRes.m8 tmp

To perform searches using GPU acceleration, you can add the `--gpu` flag to `easy-search` (see [GPU-accelerated search](https://github.com/soedinglab/MMseqs2/wiki#gpu-accelerated-search)). For frequently searched databases, you can precompute a GPU-compatible database using makepaddedseqdb, which can also be used for CPU-based searches. By default, MMseqs2 uses all visible GPUs; you can adjust this behavior using the `CUDA_VISIBLE_DEVICES` environment variable.

    mmseqs createdb examples/DB.fasta targetDB
    mmseqs makepaddedseqdb targetDB targetDB_padded
    mmseqs easy-search examples/QUERY.fasta targetDB_padded alnRes.m8 tmp --gpu 1
        
The `databases` workflow provides download and setup procedures for many public reference databases, such as the Uniref, NR, NT, PFAM and many more (see [Downloading databases](https://github.com/soedinglab/mmseqs2/wiki#downloading-databases)). For example, to download and search against a database containing the Swiss-Prot reference proteins run: 

    mmseqs databases UniProtKB/Swiss-Prot swissprot tmp
    mmseqs easy-search examples/QUERY.fasta swissprot alnRes.m8 tmp
        
The speed and sensitivity of the `search` can be adjusted with `-s` parameter and should be adapted based on your use case (see [setting sensitivity -s parameter](https://github.com/soedinglab/mmseqs2/wiki#set-sensitivity--s-parameter)). A very fast search would use a sensitivity of `-s 1.0`, while a very sensitive search would use a sensitivity of up to `-s 7.0`. A detailed guide how to speed up searches is [here](https://github.com/soedinglab/MMseqs2/wiki#how-to-control-the-speed-of-the-search).

The output can be customized with the `--format-output` option e.g. `--format-output "query,target,qaln,taln"` returns the query and target accession and the pairwise alignments in tab separated format. You can choose many different [output columns](https://github.com/soedinglab/mmseqs2/wiki#custom-alignment-format-with-convertalis).

:exclamation: `easy-search` in default computes the sequence identity by dividing the number of identical residues by the alignment length (`numIdentical/alnLen`). However, `search` [estimates](https://github.com/soedinglab/MMseqs2/wiki#how-does-mmseqs2-compute-the-sequence-identity) the identity in default. To output real sequence identity use `--alignment-mode 3` or `-a`.

### Taxonomy
The `easy-taxonomy` workflow can be used to assign sequences taxonomical labels. It performs a search against a sequence database with taxonomy information (seqTaxDb), chooses the most representative sets of aligned target sequences according to different strategies (according to `--lca-mode`) and computes the lowest common ancestor among those.

    mmseqs createdb examples/DB.fasta targetDB
    mmseqs createtaxdb targetDB tmp
    mmseqs createindex targetDB tmp
    mmseqs easy-taxonomy examples/QUERY.fasta targetDB alnRes tmp

By default, `createtaxdb` assigns a Uniprot accession to a taxonomical identifier to every sequence and downloads the NCBI taxonomy. We also support [BLAST](https://github.com/soedinglab/MMseqs2/wiki#create-a-sequence-database-with-taxonomic-information-from-an-existing-blast-database), [SILVA](https://github.com/soedinglab/MMseqs2/wiki#create-a-sequence-database-with-taxonomic-information-for-silva) or [custom taxonomical](https://github.com/soedinglab/MMseqs2/wiki#manually-annotate-a-sequence-database-with-taxonomic-information) databases. Many common taxonomic reference databases can be easily downloaded and set up by the [`databases` workflow](https://github.com/soedinglab/mmseqs2/wiki#downloading-databases).

Read more about the [taxonomy format](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-format) and the [classification](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-assignment-using-mmseqs-taxonomy) in our user guide.

### Supported search modes

MMseqs2 provides many additional search modes:
 * Iterative sequences-profile searches (like PSI-BLAST) with the `--num-iterations` parameter
 * [Translated searches](https://github.com/soedinglab/MMseqs2/wiki#translated-sequence-searching) of nucleotides against proteins (blastx), proteins against nucleotides (tblastn) or nucleotide against nucleotide (tblastx)
 * [Iterative increasing sensitivity searches](https://github.com/soedinglab/MMseqs2/wiki#how-to-find-the-best-hit-the-fastest-way) to find only the best hits faster
 * [Taxonomic assignment](https://github.com/soedinglab/MMseqs2/wiki#taxonomy-assignment-using-mmseqs-taxonomy) using 2bLCA or LCA
 * Fast ungapped alignment searches to find [very similar sequence matches](https://github.com/soedinglab/MMseqs2/wiki#mapping-very-similar-sequences-using-mmseqs-map)
 * Very fast and sensitive searches against [profile databases such as the PFAM](https://github.com/soedinglab/MMseqs2/wiki#how-to-create-a-target-profile-database-from-pfam)
 * [Reciprocal best hits search](https://github.com/soedinglab/MMseqs2/wiki#reciprocal-best-hit-using-mmseqs-rbh)
 * [Web search API and user interface](https://github.com/soedinglab/MMseqs2-App)

Many modes can also be combined. You can, for example, do a translated nucleotide against protein profile search.

### Memory requirements
MMseqs2 minimum memory requirements for `cluster` or `linclust` is 1 byte per sequence residue, `search` needs 1 byte per target residue. Sequence databases can be compressed using the `--compress` flag, DNA sequences can be reduced by a factor of `~3.5` and proteins by `~1.7`.
   
MMseqs2 checks the available system memory and automatically divides the target database in parts that fit into memory. Splitting the database will increase the runtime slightly. It is possible to control the memory usage using `--split-memory-limit`.

### How to run MMseqs2 on multiple servers using MPI
MMseqs2 can run on multiple cores and servers using OpenMP and Message Passing Interface (MPI).
MPI assigns database splits to each compute node, which are then computed with multiple cores (OpenMP).

Make sure that MMseqs2 was compiled with MPI by using the `-DHAVE_MPI=1` flag (`cmake -DHAVE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..`). Our precompiled static version of MMseqs2 cannot use MPI. The version string of MMseqs2 will have a `-MPI` suffix, if it was built successfully with MPI support.

To search with multiple servers, call the `search` or `cluster` workflow with the MPI command exported in the RUNNER environment variable. The databases and temporary folder have to be shared between all nodes (e.g. through NFS):

    RUNNER="mpirun -pernode -np 42" mmseqs search queryDB targetDB resultDB tmp

#### Multi-node `linclust` on Slurm

The default Linclust v2 workflow runs under one persistent MPI launch. Rank 0 executes
serial workflow stages, while all ranks execute `kmermatcher` and `align2clust` in-process
using the same communicator. Do not set `RUNNER` or pass `--mpi-runner` to Linclust v2.

Requirements and recommended topology:

* Build with `-DHAVE_MPI=1` (version string gets a `-MPI` suffix). A non-MPI binary
  launched under a multi-task runner (e.g. `srun -n 4 mmseqs ...` without `-DHAVE_MPI=1`)
  fails fast with an explicit error instead of silently letting every task race on the
  same output files.
* The input/output databases and the `--tmp-folder` must be on a shared POSIX filesystem
  visible to every node (e.g. NFS, Lustre, GPFS).
* Recommended: one MPI rank per node, with `--threads` set to the number of CPUs allocated
  to that node (OpenMP parallelizes within each rank). Running multiple ranks per node
  works but duplicates the memory-mapped sequence/prefilter DB cache footprint per rank, so
  it is not the default recommendation.

Example `sbatch` script requesting four nodes:

```sh
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

srun --mpi=pmix -n "$SLURM_NTASKS" --ntasks-per-node=1 \
    mmseqs linclust queryDB resultDB /shared/tmp --threads "$SLURM_CPUS_PER_TASK"
```

The equivalent direct launch is
`mpirun -np 4 mmseqs linclust queryDB resultDB /shared/tmp`. The exact
`srun --mpi=...` plugin name is site-specific.

Checkpoint/resume across topology changes: `align2clust` splits its candidate-generation
work into fixed-size, topology-independent chunks (tunable via
`MMSEQS_ALIGN2CLUST_CHUNK_SIZE`, default 50000) and records completed chunks under a
`_chunks` directory next to the output. If a run is interrupted, rerunning the identical
`mmseqs linclust`/`align2clust` command resumes by regenerating only the chunks that were
not yet completed and reusing all finished ones -- this works even if the rank count
changes between the interrupted and resumed run (e.g. starting on 8 nodes and resuming on
1), because chunk boundaries and ownership only depend on the (fixed) input size and chunk
size, never on the runner or rank count. The `_chunks` directory is removed automatically
once a run completes successfully; keep it (e.g. via `--remove-tmp-files 0`) if you want to
inspect or reuse it after a failure.

Known limitations and scaling expectations:

* `--include-align-files` is currently not supported together with distributed (multi-rank)
  candidate generation in `align2clust` and is rejected at startup in that combination; it
  still works normally on a single rank.
* Linclust v1 does not support a multi-rank outer launch. It retains the legacy
  single-rank workflow.
* `clusthash`, `clust`, and the final merge/representative-switching stages are not
  MPI-aware and run only on rank 0 while the other ranks wait.
* Representative selection and cluster assignment are a single deterministic, ordered
  reduction performed on rank 0 after all candidate chunks are collected. Scaling is
  therefore strongest when candidate/alignment generation dominates total runtime; very
  large numbers of nodes will eventually be bottlenecked by this final reduction step and by
  shared-filesystem bandwidth for chunk I/O rather than by compute.

## Contributors

MMseqs2 exists thanks to all the people who contribute. 
<a href="https://github.com/soedinglab/mmseqs2/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=soedinglab/mmseqs2" />
</a>
