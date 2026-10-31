configtuples.hpp defines kernel tile size parameters. Kernels will be compiled for the cartesian product of all parameter sets


tuning_gapless prints performance of all measured kernels in CSV format to stdout.

`gpuArch,alphabetSize,blocksize,tilesize,groupsize,numItems,mingcups,avggcups,maxgcups,datatype,approach,resultkind,numTiles`


Parameters for tuning_gapless:

    --int8: Run uint8x4 kernels

    --short2: Run s16x2 kernels

    --int32: Run int32 kernels

    --half2: Run fp16x2 kernels

    --float: Run fp32 kernels

    --endpos: Run kernels that compute both score and endpositions instead of score only

    --printBestConfigs: In addition to default output, print the best measured kernel for each tilesize in CSV format, 
                       and the corresponding C++ tuple of configs.
                       
    --tiles N: For a given tile size derived from the kernel config, use a query length of N * tilesize. subject length is always set to tilesize (default 1)

    --timingIterations N: Measure each kernel N times (default 10)

    --sleep N: Sleep N milliseconds after each config (default 100)

    --checkResults: Compare GPU alignment results with recomputed CPU results.