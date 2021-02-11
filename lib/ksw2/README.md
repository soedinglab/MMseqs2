## Introduction

KSW2 is a library to align a pair of biological sequences based on dynamic
programming (DP). So far it comes with global alignment and alignment extension
(no local alignment yet) under an affine gap cost function: gapCost(*k*) =
*q*+*k*\**e*, or a two-piece affine gap cost: gapCost2(*k*) = min{*q*+*k*\**e*,
*q2*+*k*\**e2*}. For the latter cost function, if *q*+*e*<*q2*+*e2* and *e*>*e2*,
(*q*,*e*) is effectively applied to short gaps only, while (*q2*,*e2*) applied
to gaps no shorter than ceil((*q2*-*q*)/(*e*-*e2*)-1). It helps to retain long
gaps. The algorithm behind the two-piece cost is close to [Gotoh
(1990)][piece-affine].

KSW2 supports fixed banding and optionally produces alignment paths (i.e.
CIGARs) with gaps either left- or right-aligned. It provides implementations
using SSE2 and SSE4.1 intrinsics based on [Hajime Suzuki][hs]'s diagonal
[formulation][hs-eq] which enables 16-way SSE parallelization for the most part
of the inner loop, regardless of the maximum score of the alignment.

KSW2 implements the Suzuki-Kasahara algorithm and is a component of
[minimap2][mm2]. If you use KSW2 in your work, please cite:

> * Suzuki, H. and Kasahara, M. (2018). Introducing difference recurrence relations for faster semi-global alignment of long sequences. *BMC Bioinformatics*, **19**:45.
> * Li, H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, **34**:3094-3100.

## Usage

Each `ksw2_*.c` file implements a single function and is independent of each
other. Here are brief descriptions about what each file implements:

* [ksw2_gg.c](ksw2_gg.c): global alignment; Green's standard formulation
* [ksw2_gg2.c](ksw2_gg2.c): global alignment; Suzuki's diagonal formulation
* [ksw2_gg2_sse.c](ksw2_gg2_sse.c): global alignment with SSE intrinsics; Suzuki's
* [ksw2_extz.c](ksw2_extz.c): global and extension alignment; Green's formulation
* [ksw2_extz2_sse.c](ksw2_extz2_sse.c): global and extension with SSE intrinsics; Suzuki's
* [ksw2_extd.c](ksw2_extd.c): global and extension alignment, dual gap cost; Green's formulation
* [ksw2_extd2_sse.c](ksw2_extd2_sse.c): global and extension, dual gap cost, with SSE intrinsics; Suzuki's

Users are encouraged to copy the header file `ksw2.h` and relevant
`ksw2_*.c` file to their own source code trees. On x86 CPUs with SSE2
intrinsics, `ksw2_extz2_sse.c` is recommended in general. It supports global
alignment, alignment extension with Z-drop, score-only alignment, global-only
alignment and right-aligned CIGARs. `ksw2_gg*.c` are mostly for demonstration
and comparison purposes. They are annotated with more comments and easier to
understand than `ksw2_ext*.c`. Header file [ksw2.h](ksw2.h) contains brief
documentations. TeX file [ksw2.tex](tex/ksw2.tex) gives brief derivation.

To compile the test program `ksw-test`, just type `make`. It takes the
advantage of SSE4.1 when available. To compile with SSE2 only, use `make
sse2=1` instead. If you have installed [parasail][para], use `make
parasail=prefix`, where `prefix` points to the parasail install directory (e.g.
`/usr/local`).

The following shows a complete example about how to use the library.
```c
#include <string.h>
#include <stdio.h>
#include "ksw2.h"

void align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(tseq), ql = strlen(qseq);
	uint8_t *ts, *qs, c[256];
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(c, 4, 256);
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
	for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
		printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
	putchar('\n');
	free(ez.cigar); free(ts); free(qs);
}

int main(int argc, char *argv[])
{
	align("ATAGCTAGCTAGCAT", "AGCTAcCGCAT", 1, -2, 2, 1);
	return 0;
}
```

## Performance Analysis

The following table shows timing on two pairs of long sequences (both in the
"test" directory).

|Data set|Command line options             |Time (s)|CIGAR|Ext|SIMD|Source  |
|:-------|:--------------------------------|:-------|:---:|:-:|:--:|:-------|
|50k     |-t gg -s                         |7.3     |N    |N  |N   |ksw2    |
|        |-t gg2 -s                        |19.8    |N    |N  |N   |ksw2    |
|        |-t extz -s                       |9.2     |N    |Y  |N   |ksw2    |
|        |-t ps\_nw                        |9.8     |N    |N  |N   |parasail|
|        |-t ps\_nw\_striped\_sse2\_128\_32|2.9     |N    |N  |SSE2|parasail|
|        |-t ps\_nw\_striped\_32           |2.2     |N    |N  |SSE4|parasail|
|        |-t ps\_nw\_diag\_32              |3.0     |N    |N  |SSE4|parasail|
|        |-t ps\_nw\_scan\_32              |3.0     |N    |N  |SSE4|parasail|
|        |-t extz2\_sse -sg                |0.96    |N    |N  |SSE2|ksw2    |
|        |-t extz2\_sse -sg                |0.84    |N    |N  |SSE4|ksw2    |
|        |-t extz2\_sse -s                 |3.0     |N    |Y  |SSE2|ksw2    |
|        |-t extz2\_sse -s                 |2.7     |N    |Y  |SSE4|ksw2    |
|16.5k   |-t gg -s                         |0.84    |N    |N  |N   |ksw2    |
|        |-t gg                            |1.6     |Y    |N  |N   |ksw2    |
|        |-t gg2                           |3.3     |Y    |N  |N   |ksw2    |
|        |-t extz                          |2.0     |Y    |Y  |N   |ksw2    |
|        |-t extz2\_sse                    |0.40    |Y    |Y  |SSE4|ksw2    |
|        |-t extz2\_sse -g                 |0.18    |Y    |N  |SSE4|ksw2    |

The standard DP formulation is about twice as fast as Suzuki's diagonal
formulation (`-tgg` vs `-tgg2`), but SSE-based diagonal formulation
is several times faster than the standard DP. If we only want to compute one
global alignment score, we can use 16-way parallelization in the entire inner
loop.  For extension alignment, though, we need to keep an array of 32-bit
scores and have to use 4-way parallelization for part of the inner loop. This
significantly reduces performance (`-sg` vs `-s`).  KSW2 is faster than
parasail partly because the former uses one score for all matches and another
score for all mismatches. For diagonal formulations, vectorization is more
complex given a generic scoring matrix.

It is possible to further accelerate global alignment with dynamic banding as
is implemented in [edlib][edlib]. However, it is not as effective for extension
alignment. Another idea is [adaptive banding][adap-band], which might be worth
trying at some point.

## Alternative Libraries

|Library         |CIGAR|Intra-seq|Affine-gap|Local    |Global   |Glocal   |Extension|
|:---------------|:---:|:-------:|:--------:|:-------:|:-------:|:-------:|:-------:|
|[edlib][edlib]  |Yes  |Yes      |No        |Very fast|Very fast|Very fast|N/A      |
|[KSW][klib]     |Yes  |Yes      |Yes       |Fast     |Slow     |N/A      |Slow     |
|KSW2            |Yes  |Yes      |Yes/dual  |N/A      |Fast     |N/A      |Fast     |
|[libgaba][gaba] |Yes  |Yes      |Yes       |N/A?     |N/A?     |N/A?     |Fast     |
|[libssa][ssa]   |No   |No?      |Yes       |Fast     |Fast     |N/A      |N/A      |
|[Opal][opal]    |No   |No       |Yes       |Fast     |Fast     |Fast     |N/A      |
|[Parasail][para]|No   |Yes      |Yes       |Fast     |Fast     |Fast     |N/A      |
|[SeqAn][seqan]  |Yes  |Yes      |Yes       |Slow     |Slow     |Slow     |N/A      |
|[SSW][ssw]      |Yes  |Yes      |Yes       |Fast     |N/A      |N/A      |N/A      |
|[SWIPE][swipe]  |Yes  |No       |Yes       |Fast     |N/A?     |N/A?     |N/A      |
|[SWPS3][swps3]  |No   |Yes      |Yes       |Fast     |N/A?     |N/A      |N/A      |



[hs]: https://github.com/ocxtal
[hs-eq]: https://github.com/ocxtal/diffbench
[edlib]: https://github.com/Martinsos/edlib
[klib]: https://github.com/attractivechaos/klib
[para]: https://github.com/jeffdaily/parasail
[opal]: https://github.com/Martinsos/opal
[ssw]: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
[ssa]: https://github.com/RonnySoak/libssa
[gaba]: https://github.com/ocxtal/libgaba
[adap-band]: https://github.com/ocxtal/adaptivebandbench
[swipe]: https://github.com/torognes/swipe
[swps3]: http://lab.dessimoz.org/swps3/
[seqan]: http://seqan.de
[piece-affine]: https://www.ncbi.nlm.nih.gov/pubmed/2165832
[mm2]: https://github.com/lh3/minimap2
