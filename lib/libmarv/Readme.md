




# CUDASW++4.0-GapLessFilter


## Software requirements
* Linux operating system with compatible CUDA Toolkit 12 or newer
* C++17 compiler
* zlib
* make

## Hardware requirements
*   A modern CUDA-capable GPU of generation Ampère or newer. We have tested CUDASW4 on Ampère (sm_80), Ada Lovelace (sm_89), and Hopper (sm_90). Older generations lack hardware-support for specific instructions and may run at reduced speeds or may not run at all.




## Build
Our software has two components, **makedb** and **align** . **makedb** is used to construct a database which can be queried by **align**.

The build step compiles the GPU code for all GPU archictectures of GPUs detected in the system. The CUDA environment variable `CUDA_VISIBLE_DEVICES` can be used to control the detected GPUs. If `CUDA_VISIBLE_DEVICES` is not set, it will default to all GPUs in the system.

* Build makedb: `make makedb`

* Build align: `make align`

* Build align for the GPU architecture of GPUs 0 and 1: `CUDA_VISIBLE_DEVICES=0,1 make align`

## Database construction
Use **makedb** to create a database from a fasta file. The file can be gzip'ed.
We support fasta files with up to 2 billion sequences.

```
mkdir -p dbfolder
./makedb input.fa(.gz) dbfolder/dbname [options]
```

Options:
* --mem val : Memory limit. Can use suffix K,M,G. If makedb requires more memory, temp files in temp directory will be used. Default all available memory.
* --tempdir val : Temp directory for temporary files. Must exist. Default is db output directory.



## Querying the database
Use **align** to query the database. **align** has two mandatory arguments. 
1. `--query` The query file which contains all queries
2. `--db` The path to the reference database constructed with makedb. 

Run `./align --help` to get a complete list of options.

By default, the results will be output to stdout in plain text. Results can be output to file instead (`--of filename`), and can be output as tab-separated values (`--tsv`). Example tsv output is given below.

| Query number | Query length | Query header | Result number | Result score | Reference length | Reference header | Reference ID in DB | Alignment_end_query | Alignment_end_ref |
|------------|------------|------------|------------|------------|------------| ------------|------------|------------|------------|
| 0 | 144 | gi\|122087146 | 0 | 541 | 148 | UniRef50_P02233 | 23128215 | 100 | 100 |
| 0 | 144 | gi\|122087146 | 1 | 444 | 144 | UniRef50_P02238  | 22381647 | 100 | 100 |


## Selecting GPUs
Similar to the build process, **align** will use all GPUs that are set with `CUDA_VISIBLE_DEVICES`, or all GPUs if `CUDA_VISIBLE_DEVICES` is not set. 

```
# use the gpus that are currently set in CUDA_VISIBLE_DEVICES
./align --query queries.fa(.gz) --db dbfolder/dbname

# use gpus 0 and 1 for only this command
CUDA_VISIBLE_DEVICES=0,1 ./align --query queries.fa(.gz) --db dbfolder/dbname
```

## Scoring options

```
    --top val : Output the val best scores. Default val = 10.
    --gop val : Gap open score. Overwrites blosum-dependent default score.
    --gex val : Gap extend score. Overwrites blosum-dependent default score.
    --mat val : Set substitution matrix. Supported values: blosum45, blosum50, blosum62, blosum80. Default val = blosum62.
    --scanType val : Set scan type. Supported values = {Gapless, SW_Endpos, Gapless+SW_Endpos}.
            Gapless: Scan whole DB with gapless alignment. 
            SW_Endpos: Scan whole DB with Smith Waterman Alignment, output score and end position.
            Gapless+SW_Endpos: Scan whole DB with gapless alignment, then re-scan top results with Smith Waterman. Default val = Gapless
    --subjectIdsFile val : Only consider database sequences with index specified in file. Must be a text file, one index per line.
            Do not use together with scanType Gapless+SW_Endpos. When --subjectIdsFile is set, option --top is ignored.
```


## Memory options

```
    --maxGpuMem val : Try not to use more than val bytes of gpu memory per gpu. This is not a hard limit. Can use suffix K,M,G. All available gpu memory by default.
    --maxTempBytes val : Size of temp storage in GPU memory. Can use suffix K,M,G. Default val = 4G
    --maxBatchBytes val : Process DB in batches of at most val bytes. Can use suffix K,M,G. Default val = 128M
    --maxBatchSequences val : Process DB in batches of at most val sequences. Default val = 10000000
```

Depending on the database size and available total GPU memory, the database is transferred to the GPU once for all queries, or it is processed in batches which requires a transfer for each query. Above options give some control over memory usage. For best performance, the complete database must fit into `maxGpuMem` times the number of used GPUs.

## Other options
```
    --verbose : More console output. Shows timings.
    --printLengthPartitions : Print number of sequences per length partition in db.
    --interactive : Loads DB, then waits for sequence input by user
    --help : Print all options
```

## Peak performance test:

Align all queries to a simulated database of 2000000 sequences of length 512. All sequences are identical.
```
./align --query allqueries.fasta --top 0 --verbose --uploadFull --pseudodb 2000000 512 1
```

Align all queries to a simulated database of 2000000 sequences of length 512. All sequences are random.
```
./align --query allqueries.fasta --top 0 --verbose --uploadFull --pseudodb 2000000 512 0
```