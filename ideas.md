# Old ideas and history
Originally, this project started with trying to improve Daily's prefix scan algorithm.
The goal was to make it (static) banded and also use a difference encoding like in the
Suzuki-Kasahara algorithm. Additionally, it would be nice to support protein alignment (in addition to
DNA alignment) but with very narrow lanes, for maximum SIMD parallelism.

This plan did not really work out, so I looked into pivoting to adaptive banding methods.
Adaptive banding allows a very small band to be used, compared to the traditional static
banding approach. Part of this was to stay competitive with the recent Wavefront Alignment
algorithm in terms of speed, and also improve on it since it could not handle complex amino
acid scoring schemes that was necessary in protein alignment.

However, with vertical or horizontal bands, I quickly realized that it was too hard to
identify the direction to shift. I also thought about using L-shaped areas and other shapes
to tile the DP matrix, but eventually I settled on square blocks due to their flexibility.

WASM SIMD support was especially interesting since it is cross-platform and runs in the
browser. Who knows? Maybe people will start developing more bioinformatics tools that
run in the browser!

Compared to other algorithms, block aligner is very optimistic. Try small blocks before larger blocks!

## Some failed ideas
1. What if we took Daily's prefix scan idea and made it faster and made it banded using
ring buffers and had tons of 32-bit offsets for intervals of the band to prevent overflow?
(This actually works, but it is soooooo complex.)
2. What if we took that banded idea (a single thin vertical band) and made it adaptive?
3. What if we placed blocks like Minecraft, where there is no overlap between blocks?
4. What if we compared the rightmost column and bottommost row in each block to decide
which direction to shift? (Surprisingly, using the first couple of values in each column
or row is better than using the whole column/row. Also, comparing the sum of scores worked
better than comparing the max. Update: revisiting this, it seems like max actually works best!)
5. Use a branch-predictor-like scheme to predict which direction to shift as a tie-breaker
when shifting right or down seem equally good.
6. ...

Some old code is located in the `src/old` directory.
