//
// Created by mad on 2/3/16.
//

#include <Debug.h>
#include <Util.h>
#include "MsaFilter.h"
#include "simd.h"
#include "MathUtil.h"
#include "MultipleAlignment.h"

MsaFilter::MsaFilter(int maxSeqLen, int maxSetSize, SubstitutionMatrix *m, int gapOpen, int gapExtend) :
    // TODO allow changing these?
    PLTY_GAPOPEN(6.0f), PLTY_GAPEXTD(1.0f), gapOpen(gapOpen), gapExtend(gapExtend) {
    this->m = m;
    this->maxSeqLen = maxSeqLen;
    this->maxSetSize = maxSetSize;
    this->Nmax = new int[maxSeqLen + 2];       // position-dependent maximum-sequence-identity threshold for filtering? (variable used in former version was idmax)
    this->idmaxwin = new int[maxSeqLen + 2];   // minimum value of idmax[i-WFIL,i+WFIL]
    this->N = new int[maxSeqLen + 2];  // N[i] number of already accepted sequences at position i
    this->in = (char*)malloc((maxSetSize + 1) * sizeof(char));  // in[k]=1: seq k has been accepted; in[k]=0: seq k has not yet been accepted at current seqid
    this->inkk = (char*)malloc((maxSetSize + 1) * sizeof(char));  // inkk[k]=1 iff in[ksort[k]]=1 else 0;
    this->seqid_prev = (int*)malloc((maxSetSize + 1) * sizeof(int));  // maximum-sequence-identity threshold used in previous round of filtering (with lower seqid)
    this->first = (int*)malloc(maxSetSize * sizeof(int));         // first non-gap position in sequence k
    this->last = (int*)malloc(maxSetSize * sizeof(int));          // last  non-gap position in sequence k
    this->nres = (int*)malloc(maxSetSize * sizeof(int));
    this->ksort = (int*)malloc(maxSetSize * sizeof(int));
    this->display = (char*)malloc((maxSetSize + 2) * sizeof(char));
    this->keep = (char*)malloc(maxSetSize * sizeof(char));
}

MsaFilter::~MsaFilter() {
    delete [] Nmax;
    delete [] idmaxwin;
    delete [] N;
    free(in);
    free(inkk);
    free(seqid_prev);
    free(first);
    free(last);
    free(nres);
    free(ksort);
    free(display);
    free(keep);
}

void MsaFilter::increaseSetSize(int newSetSize) {
    if (newSetSize > maxSetSize) {
        maxSetSize = newSetSize * 1.5;
        in = (char*)realloc(in, maxSetSize * sizeof(char));
        inkk = (char*)realloc(inkk, maxSetSize * sizeof(char));
        seqid_prev = (int*)realloc(seqid_prev, maxSetSize * sizeof(int));
        first = (int*)realloc(first, maxSetSize * sizeof(int));
        last = (int*)realloc(last, maxSetSize * sizeof(int));
        nres = (int*)realloc(nres, maxSetSize * sizeof(int));
        ksort = (int*)realloc(ksort, maxSetSize * sizeof(int));
        display = (char*)realloc(display, maxSetSize * sizeof(char));
        keep = (char*)realloc(keep, maxSetSize * sizeof(char));
    }
}

size_t MsaFilter::filter(MultipleAlignment::MSAResult &msa, std::vector<Matcher::result_t> &alnResults, int coverage, int qid, float qsc, int max_seqid, int Ndiff) {
    size_t filteredSize = filter(msa.setSize, msa.centerLength, coverage, qid, qsc, max_seqid, Ndiff, (const char **) msa.msaSequence, true);
    if (!alnResults.empty()) {
        // alignmentResults does not include the query
        for (size_t i = 0, j = 0; j < msa.setSize - 1; j++) {
            if (keep[j + 1] != 0) {
                if (i < j) {
                    std::swap(alnResults[i], alnResults[j]);
                }
                i++;
            }
        }
        alnResults.resize(filteredSize - 1);
    }
    return filteredSize;
}

size_t MsaFilter::filter(const int N_in, const int L, const int coverage, const int qid,
                       const float qsc, const int max_seqid, int Ndiff, const char **X, const bool shuffleMsa) {
    increaseSetSize(N_in);

    int seqid1 = 20;
    // X[k][i] contains column i of sequence k in alignment (first seq=0, first char=1) (0-3: ARND ..., 20:X, 21:GAP)
//    char** X = (char **) &msaSequence;

    // In the beginnning, keep[k] is 1 for all regular amino acid sequences and 0 for all others (ss_conf, ss_pred,...)
    // In the end, keep[k] will be 1 for all regular representative sequences kept in the alignment, 0 for all others
    // Sequences with keep[k] = 2 will cannot be filtered out and will remain in the alignment.
    // If a consensus sequence exists it has k = kfirst and keep[k] = 0, since it should not enter into the profile calculation.
    const int WFIL = 25;            // see previous line

    int diffNmax = Ndiff;    // current  maximum difference of Nmax[i] and Ndiff
    int diffNmax_prev = 0;   // previous maximum difference of Nmax[i] and Ndiff
    int seqid;  // current  maximum value for the position-dependent maximum-sequence-identity thresholds in idmax[]
    int seqid_step = 0;         // previous increment of seqid

    float diff_min_frac;  // minimum fraction of differing positions between sequence j and k needed to accept sequence k
    float qdiff_max_frac = 0.9999 - 0.01 * qid;  // maximum allowable number of residues different from query sequence
    int diff = 0;  // number of differing positions between sequences j and k (counted so far)
    int diff_suff;  // number of differing positions between sequences j and k that would be sufficient
    int qdiff_max;  // maximum number of residues required to be different from query
    int cov_kj;  // upper limit of number of positions where both sequence k and j have a residue
    int first_kj;             // first non-gap position in sequence j AND k
    int last_kj;              // last  non-gap position in sequence j AND k
    int kk, jj;               // indices for sequence from 1 to N_in
    int k, j;                 // kk=ksort[k], jj=ksort[j]
    int i;                    // counts residues
    int n;                    // number of sequences accepted so far
    int kfirst = 0;           // index of first real sequence

    // map data to X
    for (k = 0; k < N_in; ++k) {
        // sequence 0 is the center (query)
        keep[k] = (k == 0) ? 2 : 1;
    }
    // Initialize in[k]
    for (n = k = 0; k < N_in; ++k) {
        if (keep[k] == 2) {
            in[k] = 2;
            n++;
        } else {
            in[k] = 0;
        }
    }
    // Determine first[k], last[k]?
    for (k = 0; k < N_in; ++k)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
    {
        for (i = 0; i < L; ++i)
            if (X[k][i] < MultipleAlignment::NAA)
                break;
        first[k] = i;
        for (i = (L - 1); i > 0; i--)
            if (X[k][i] < MultipleAlignment::NAA)
                break;
        last[k] = i;
    }

    // Determine number of residues nres[k]?
    for (k = 0; k < N_in; ++k)  // do this for ALL sequences, not only those with in[k]==1 (since in[k] may be display[k])
    {
        int nr = 0;
        for (i = first[k]; i <= last[k]; ++i)
            if (X[k][i] < MultipleAlignment::NAA)
                nr++;
        this->nres[k] = nr;
//        printf("%d nres=%3i  first=%3i  last=%3i\n",k,nr,first[k],last[k]);
        if (nr == 0)
            keep[k] = 0;
    }

    std::pair<int, int> * tmpSort = new std::pair<int, int>[N_in];
    // create sorted index according to length (needed for the pairwise seq. id. comparision); afterwards, nres[ksort[kk]] is sorted by size
    for (k = 0; k < N_in; ++k) {
        tmpSort[k].first = nres[k];
        tmpSort[k].second = k;
    }
    //Sort sequences after query (first sequence) in descending order
    struct sortPairDesc {
        bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) const {
            return left.first > right.first;
        }
    };
    std::stable_sort(tmpSort + 1, tmpSort + N_in, sortPairDesc());
    for (k = 0; k < N_in; ++k) {
        ksort[k] =  tmpSort[k].second;
    }
    delete [] tmpSort;

    for (kk = 0; kk < N_in; ++kk) {
        inkk[kk] = in[ksort[kk]];
    }

    // Initialize N[i], idmax[i], idprev[i]
    for (i = 0; i < first[kfirst]; ++i)
        N[i] = 0;
    for (i = first[kfirst]; i <= last[kfirst]; ++i)
        N[i] = 1;
    for (i = last[kfirst] + 1; i < L; ++i)
        N[i] = 0;
    for (i = 0; i < L; ++i) {
        Nmax[i] = 0;
        idmaxwin[i] = -1;
    }
    for (k = 0; k < N_in; ++k)
        seqid_prev[k] = -1;
    if (Ndiff <= 0 || Ndiff >= N_in) {
        seqid1 = max_seqid;
        Ndiff = N_in;
        diffNmax = Ndiff;
    }

    // Check coverage and sim-to-query criteria for each sequence k
    for (k = 0; k < N_in; ++k) {
        if (keep[k] == 0 || keep[k] == 2)
            continue;  // seq k not regular sequence OR is marked sequence
        if (100 * nres[k] < coverage * L) {
            keep[k] = 0;
            continue;
        }  // coverage too low? => reject once and for all

        float qsc_sum = 0.0;

        // Check if score-per-column with query is at least qsc
        if (qsc > -10) {
            float qsc_min = qsc * nres[k];  // minimum total score of seq k with query

            int gapq = 0, gapk = 0;  // number of consecutive gaps in query or k'th sequence at position i
            for (int i = first[k]; i <= last[k]; ++i) {
                if (X[k][i] < 20) {
                    gapk = 0;
                    if (X[kfirst][i] < 20) {
                        gapq = 0;
                        qsc_sum += static_cast<float>(m->subMatrix[(int) X[kfirst][i]][(int) X[k][i]]);
                    } else if (X[kfirst][i] == MultipleAlignment::ANY)
                        // Treat score of X with other amino acid as 0.0
                        continue;
                    else if (gapq++)
                        qsc_sum -= PLTY_GAPEXTD;
                    else
                        qsc_sum -= PLTY_GAPOPEN;
                } else if (X[k][i] == MultipleAlignment::ANY)
                    // Treat score of X with other amino acid as 0.0
                    continue;
                else if (X[kfirst][i] < 20) {
                    gapq = 0;
                    if (gapk++)
                        qsc_sum -= PLTY_GAPEXTD;
                    else
                        qsc_sum -= PLTY_GAPOPEN;
                }
            }
//        printf("k=%3i qsc=%6.2f\n",k,qsc_sum);
            if (qsc_sum < qsc_min) {
                keep[k] = 0;
                continue;
            }  // too different from query? => reject once and for all
        }
//                  printf("  diff=%4i\n",diff);
        //Check if sequence similarity with query at least qid?
        if (qdiff_max_frac < 0.999) {

            qdiff_max = int(qdiff_max_frac * nres[k] + 0.9999);
//                  printf("k=%-4i  nres=%-4i  qdiff_max=%-4i first=%-4i last=%-4i",k,nres[k],qdiff_max,first[k],last[k]);
            diff = 0;
            for (int i = first[k]; i <= last[k]; ++i)
                // enough different residues to reject based on minimum qid with query? => break
                if (X[k][i] < MultipleAlignment::NAA
                    && X[k][i] != X[kfirst][i] && ++diff >= qdiff_max)
                    break;
//                  printf("  diff=%4i\n",diff);
            if (diff >= qdiff_max) {
                keep[k] = 0;
                continue;
            }  // too different from query? => reject once and for all
        }
        // printf("%d  qsc=%6.2f  %d  %d qid=%6.2f  \n",k,qsc_sum/nres[k],nres[k], diff, 100.0*(1.0-(float)(diff)/nres[k]));
    }

    // If no sequence left, issue warning and put back first real sequence into alignment
    int nn=0;
    for (k = 0; k < N_in; ++k) {
        if (keep[k] > 0) {
            nn++;
        }
    }

    if (nn == 0) {
        for (k = 0; k < N_in; k++) {
            if (display[k] != 2) {
                keep[k] = 1;
                break;
            }
        }
        if (keep[k] == 1) {
//            Debug(Debug::WARNING) << "Warning in " << __FILE__ << ":" << __LINE__
//            << ": " << __func__ << ":" << "\n";
//            Debug(Debug::WARNING)
//            << "\tFiltering removed all sequences in alignment. Inserting back first sequence.\n";
            ;
        } else if (display[kfirst] == 2) {  // the only sequence in the alignment is the consensus sequence :-(
//            Debug(Debug::WARNING) << "Warning in " << __FILE__ << ":" << __LINE__
//            << ": " << __func__ << ":" << "\n";
//            Debug(Debug::WARNING)
//            << "\tAlignment contains no sequence except consensus sequence. Using consensus sequence for searching.\n";
            ;
        } else {
            Debug(Debug::WARNING) << "The alingment %s does not contain any sequences.\n";
        }
    }

    // If min required seqid larger than max required seqid, return here without doing pairwise seqid filtering
    if (seqid1 > max_seqid) {
        if (shuffleMsa) {
            shuffleSequences(X, N_in);
        }
        return nn;
    }

    // Successively increment idmax[i] at positons where N[i]<Ndiff
    seqid = seqid1;
    while (seqid <= max_seqid) {
        bool stop = true;
        // Update Nmax[i]
        diffNmax_prev = diffNmax;
        diffNmax = 0;
        for (i = 0; i < L; ++i) {
            int max = 0;
            for (j = std::max(0, std::min(L - 2 * WFIL + 1, i - WFIL));
                 j < std::min(L, std::max(2 * WFIL, i + WFIL)); ++j)
                if (N[j] > max)
                    max = N[j];
            if (Nmax[i] < max)
                Nmax[i] = max;
            if (Nmax[i] < Ndiff) {
                stop = false;
                idmaxwin[i] = seqid;
                if (diffNmax < Ndiff - Nmax[i])
                    diffNmax = Ndiff - Nmax[i];
            }
        }

//        printf("seqid=%3i  diffNmax_prev= %-4i   diffNmax= %-4i   n=%-5i  N_in-N_ss=%-5i\n",seqid,diffNmax_prev,diffNmax,n,N_in);
        if (stop) {
            break;
        }

//       // DEBUG
//       printf("idmax    ");
//       for (i=1; i<=L; ++i) printf("%2i ",idmax[i]);
//       printf("\n");
//       printf("idmaxwin ");
//       for (i=1; i<=L; ++i) printf("%2i ",idmaxwin[i]);
//       printf("\n");
//       printf("N[i]     ");
//       for (i=1; i<=L; ++i) printf("%2i ",N[i]);
//       printf("\n");

        // Loop over all candidate sequences kk (-> k)
        for (kk = 0; kk < N_in; ++kk) {
            if (inkk[kk])
                continue;   // seq k already accepted
            k = ksort[kk];
            if (!keep[k])
                continue;  // seq k is not regular aa sequence or already suppressed by coverage or qid criterion
            if (keep[k] == 2) {
                inkk[kk] = 2;
                continue;
            }  // accept all marked sequences (no n++, since this has been done already)

            // Calculate max-seq-id threshold seqidk for sequence k (as maximum over idmaxwin[i])
            if (seqid >= 100) {
                in[k] = inkk[kk] = 1;
                n++;
                continue;
            }

            float seqidk = seqid1;
            for (i = first[k]; i <= last[k]; ++i)
                if (idmaxwin[i] > seqidk)
                    seqidk = idmaxwin[i];
            if (seqid == seqid_prev[k])
                continue;  // sequence has already been rejected at this seqid threshold => reject this time
            seqid_prev[k] = seqid;
            diff_min_frac = 0.9999 - 0.01 * seqidk;  // min fraction of differing positions between sequence j and k needed to accept sequence k
            // Loop over already accepted sequences
            for (jj = 0; jj < kk; ++jj) {
                if (!inkk[jj])
                    continue;
                j = ksort[jj];
//                for (int i = first[k]; i <= last[k]; ++i) {
//                    printf("%02d", (int)X[k][i]);
//                }
//                std::cout << std::endl;
//                for (int i = first[j]; i <= last[j]; ++i) {
//                    printf("%02d", (int)X[j][i]);
//                }
//                std::cout << std::endl;

                first_kj = std::max(first[k], first[j]);
                last_kj = std::min(last[k], last[j]);
                cov_kj = last_kj - first_kj + 1;
                diff_suff = int(diff_min_frac * std::min(nres[k], cov_kj) + 0.999);  // nres[j]>nres[k] anyway because of sorting
                diff = 0;
                const simd_int * XK = (simd_int *) X[k];
                const simd_int * XJ = (simd_int *) X[j];
                const int first_kj_simd = first_kj / (VECSIZE_INT * 4);
                const int last_kj_simd = last_kj / (VECSIZE_INT * 4) + 1;
                // coverage correction for simd
                // because we do not always hit the right start with simd.
                // This works because all sequence vector are initialized with GAPs so the sequnces is surrounded by GAPs
                const int first_diff_simd_scalar = std::abs(
                        first_kj_simd * (VECSIZE_INT * 4) - first_kj);
                const int last_diff_simd_scalar = std::abs(
                        last_kj_simd * (VECSIZE_INT * 4) - (last_kj + 1));

                cov_kj += (first_diff_simd_scalar + last_diff_simd_scalar);

                // _mm_set1_epi8 pseudo-instruction is slow!
                const simd_int NAAx16 = simdi8_set(MultipleAlignment::NAA - 1);
                for (int i = first_kj_simd; i < last_kj_simd && diff < diff_suff; ++i) {
                    // None SIMD function
                    // enough different residues to accept? => break
                    // if (X[k][i] >= NAA || X[j][i] >= NAA)
                    //    cov_kj--;
                    // else if (X[k][i] != X[j][i] && ++diff >= diff_suff)
                    //    break; // accept (k,j)

                    const simd_int NO_AA_K = simdi8_gt(XK[i], NAAx16);  // pos without amino acid in seq k
                    const simd_int NO_AA_J = simdi8_gt(XJ[i], NAAx16);  // pos without amino acid in seq j

                    // Compute 16 bits indicating positions with GAP, ANY or ENDGAP in seq k or j
                    // int _mm_movemask_epi8(__m128i a) creates 16-bit mask from most significant bits of
                    // the 16 signed or unsigned 8-bit integers in a and zero-extends the upper bits.
                    int res = simdi8_movemask(simdi_or(NO_AA_K, NO_AA_J));
//                    for (int u = 0; u < 32; ++u) {
//                        printf("%02d:%02d ", (int) ((char*)&XK[i])[u], (int) ((char*)&XK[i])[u]);
//                    }
//                    std::cout << std::endl;
                    cov_kj -= __builtin_popcount(res);  // subtract positions that should not contribute to coverage

                    // Compute 16 bit mask that indicates positions where k and j have identical residues
                    int c = simdi8_movemask(simdi8_eq(XK[i], XJ[i]));

                    // Count positions where  k and j have different amino acids, which is equal to 16 minus the
                    //  number of positions for which either j and k are equal or which contain ANY, GAP, or ENDGAP
                    diff += (VECSIZE_INT * 4) - __builtin_popcount(c | res);
                }
//            // DEBUG
//            printf("%20.20s with %20.20s:  diff=%i  diff_min_frac*cov_kj=%f  diff_suff=%i  nres=%i  cov_kj=%i\n",sname[k],sname[j],diff,diff_min_frac*cov_kj,diff_suff,nres[k],cov_kj);
//            printf("%s\n%s\n\n",seq[k],seq[j]);
                if (diff < diff_suff && float(diff) <= diff_min_frac * cov_kj && cov_kj > 0)
                    break;  //dissimilarity < acceptace threshold? Reject!

            }
            if (jj >= kk)  // did loop reach end? => accept k. Otherwise reject k (the shorter of the two)
            {
                in[k] = inkk[kk] = 1;
                n++;
                for (i = first[k]; i <= last[k]; ++i)
                    N[i]++;  // update number of sequences at position i
//            printf("%i %20.20s accepted\n",k,sname[k]);
            }
//        else
//          {
//            printf("%20.20s rejected: too similar with seq %20.20s  diff=%i  diff_min_frac*cov_kj=%f  diff_suff=%i  nres=%i  cov_kj=%i\n","blub","blub",diff,diff_min_frac*cov_kj,diff_suff,nres[k],cov_kj);
//            printf("%s\n%s\n\n",msaSequence[k],msaSequence[j]);
//          }

        }  // End Loop over all candidate sequences kk

//       // DEBUG
//       printf("\n");
//       printf("seqid_prev[k]= \n");
//       for (k=0; k<N_in; ++k) printf("%2i ",seqid_prev[k]);
//       printf("\n");

        // Increment seqid
        seqid_step = std::max(
                1, std::min(5, diffNmax / (diffNmax_prev - diffNmax + 1) * seqid_step / 2));
        seqid += seqid_step;

    }  // End Loop over seqid

//    Debug(Debug::WARNING) << n << " out of " << N_in << " sequences passed filter (";
//    if (coverage) {
//        Debug(Debug::WARNING) << coverage << "% min coverage, ";
//    }
//    if (qid) {
//        Debug(Debug::WARNING) << qid << "% min sequence identity to query, ";
//    }
//    if (qsc > -10) {
//        Debug(Debug::WARNING) << qsc << " bits min score per column to query, ";
//    }
//    if (Ndiff < N_in && Ndiff > 0) {
//        Debug(Debug::WARNING)
//        << "up to " << seqid
//        << "% position-dependent max pairwise sequence identity)\n";
//    } else {
//        Debug(Debug::WARNING) << seqid1 << "% max pairwise sequence identity)\n";
//    }

    for (k = 0; k < N_in; ++k) {
        keep[k] = in[k];
    }

    if (shuffleMsa) {
        shuffleSequences(X, N_in);
    }
    return n;
}

void MsaFilter::shuffleSequences(const char ** X, size_t setSize) {
    for (size_t i = 0, j = 0; j < setSize; j++) {
        if (keep[j] != 0) {
            if (i < j) {
                const char* temp = X[i];
                X[i] = X[j];
                X[j] = temp;
            }
            i++;
        }
    }
}

void MsaFilter::getKept(bool *kept, size_t setSize) {
    for (size_t i = 0; i < setSize; i++) {
        kept[i] = keep[i] != 0;
    }
}

void MsaFilter::pruneAlignment(char ** msaSequence, int N_in, int L) {
    int   bg = 5;  // below this number of end gaps the loose HSP pruning score is used
    float bl = 0.0;   // minimum per-residue bit score with query at ends of HSP for loose end pruning
    float bs = 0.8;   // minimum per-residue bit score with query at ends of HSP for strict end pruning
    for(int seqIdx = 1; seqIdx < N_in; seqIdx++ ){
        int qfirst = 0;             // index of first query residue in pairwise alignment
        for(int i = 0; i < L; i++){
            if(msaSequence[seqIdx][i]==MultipleAlignment::GAP){
                qfirst++;
            }else{
                break;
            }
        }
        int qlast  = L - 1;         // index of last  query residue in pairwise alignment
        for(int i = qlast; i >= 0; i--){
            if(msaSequence[seqIdx][i]==MultipleAlignment::GAP){
                qlast--;
            }else{
                break;
            }
        }
        // Count gaps in template that are aligned with match residues to the left of HSP
        int gapsleft= qfirst;
        int gapsright= L - qlast;
        float bleft   = (gapsleft  >= bg) ? bs : bl;
        float bright =  (gapsright >= bg) ? bs : bl;
        int i1 = 0;
        int i2 = L;

        if( bleft > -9){
            i1 = prune(0,     L, bleft,  msaSequence[0], msaSequence[seqIdx]);
        }
        if (bright > -9)
        {
            i2 = prune(L - 1, 0, bright, msaSequence[0], msaSequence[seqIdx]);
        }
        if(i1 > 0){
            for (int i = 0; i <= i1; i++) {
                msaSequence[seqIdx][i]=MultipleAlignment::GAP;
            }
        }
        if(i2 < L - 1){
            for (int i = i2;  i< L; i++) {
                msaSequence[seqIdx][i]=MultipleAlignment::GAP;
            }
        }
    }

}


int MsaFilter::prune(int start, int end, float b, char * query, char *target) {
    float smin  = 0.0;
    float score = 0.0;
    bool gap = false;
    int i_ret = start;
    bool rev = false;
    int pos = start;
    if(end < start){
        int tmp = start;
        start = end;
        end = tmp;
        rev = true;
    }
    for(int i = start; i < end; i++) {
        if(rev == true){
            if( pos < i_ret - 20.0)
                break;
        }else{
            if( pos > i_ret + 20.0)
                break;
        }
        smin += b;
        if(query[pos] < MultipleAlignment::NAA && target[pos] < MultipleAlignment::NAA) {
            score += (static_cast<float>(m->subMatrix[(int)query[pos]][(int)target[pos]]) ) * 0.3322;
            gap = false;
        } else if(query[pos] == MultipleAlignment::GAP  || target[pos] == MultipleAlignment::GAP) {
            score -= (gap == false) ? static_cast<float>(gapOpen)* 0.3322 : static_cast<float>(gapExtend) * 0.3322;
            gap = true;
        }
        if (score < smin) {
            i_ret = pos;
            smin=0; score=0;
        }
        pos += (rev==true)? -1 : +1;
    }
    return i_ret;
}
