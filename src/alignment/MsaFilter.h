//
// Created by mad on 2/3/16.
//

#ifndef MMSEQS_MSAFILTER_H
#define MMSEQS_MSAFILTER_H


#include <SubstitutionMatrix.h>
#include "MultipleAlignment.h"

class MsaFilter {

public:

    MsaFilter(int maxSeqLen, int maxSetSize, SubstitutionMatrix *m, int gapOpen, int gapExtend);

    ~MsaFilter();
    /////////////////////////////////////////////////////////////////////////////////////
    // Select set of representative sequences in the multiple sequence alignment
    // Filter criteria:
    //   * Remove sequences with coverage of query less than "coverage" percent
    //   * Remove sequences with sequence identity to query of less than "qid" percent
    //   * If Ndiff==0, remove sequences with seq. identity larger than seqid2(=max_seqid) percent
    //   * If Ndiff>0, remove sequences with minimum-sequence-identity filter of between seqid1
    //     and seqid2 (%), where the minimum seqid threshold is determined such that,
    //     in all column blocks of at least WMIN=25 residues, at least Ndiff sequences are left.
    //     This ensures that in multi-domain proteins sequences covering one domain are not
    //     removed completely because sequences covering other domains are more diverse.
    //
    // Allways the shorter of two compared sequences is removed (=> sort sequences by length first).
    // Please note: sequence identity of sequence x with y when filtering x is calculated as
    // number of residues in sequence x that are identical to an aligned residue in y / number of residues in x
    // Example: two sequences x and y are 100% identical in their overlapping region but one overlaps by 10% of its
    // length on the left and the other by 20% on the right. Then x has 10% seq.id with y and y has 20% seq.id. with x.
    /////////////////////////////////////////////////////////////////////////////////////
    size_t filter(MultipleAlignment::MSAResult& msa, std::vector<Matcher::result_t> &alnResults, int coverage,
                  const std::vector<int> & qid_vec, float qsc, int max_seqid, int Ndiff, int filterMinEnable);
    size_t filter(const int N_in, const int L, const int coverage, const std::vector<int> & qid,
                  const float qsc, const int max_seqid, int Ndiff, int filterMinEnable, const char **X, const bool shuffleMsa);

    void getKept(bool *offsets, size_t setSize);

    const float PLTY_GAPOPEN; // for -qsc option (filter for min similarity to query): 6 bits to open gap
    const float PLTY_GAPEXTD; // for -qsc option (filter for min similarity to query): 1 bit to extend gap

    void pruneAlignment(char ** msaSequence, int N_in, int L);
	
	
private:
    // shuffles the filtered sequences to the back of the array, the unfiltered ones remain in the front
    void shuffleSequences(const char ** X, size_t setSize);

    // prune sequence based on score
    int prune(int start, int end, float b, char * query, char *target);

    void increaseSetSize(int newSetSize);

    BaseMatrix *m;

    int maxSeqLen;
    int maxSetSize;
    int gapOpen;
    int gapExtend;

    // position-dependent maximum-sequence-identity threshold for filtering? (variable used in former version was idmax)
    int *Nmax;
    // minimum value of idmax[i-WFIL,i+WFIL]
    int *idmaxwin;
    // N[i] number of already accepted sequences at position i
    int *N;
    // in[k]=1: seq k has been accepted; in[k]=0: seq k has not yet been accepted at current seqid
    char *in;
    // inkk[k]=1 iff in[ksort[k]]=1 else 0;
    char *inkk;
    // maximum-sequence-identity threshold used in previous round of filtering (with lower seqid)
    int *seqid_prev;
    int *nres;
    // first residue in sequence k
    int* first;
    // last  residue in sequence k
    int* last;
    // index for sorting sequences: X[ksort[k]]
    int* ksort;
    // display[k]=1 if sequence will be displayed in output alignments; 0 otherwise (first=0)
    char* display;
    // keep[k]=1 if sequence is included in amino acid frequencies; 0 otherwise (first=0)
    char *keep;
    // needed for qid bucket logic
    const char **X;
    char **keep_local;
};


#endif //MMSEQS_MSAFILTER_H
