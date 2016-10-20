//
// Created by mad on 10/20/16.
//

#ifndef MMSEQS_TRANSLATE_H
#define MMSEQS_TRANSLATE_H

#include <string.h>

class Translate {
public:
// translation tables specific to each genetic code instance
    char  m_AminoAcid [4097];
    char  m_OrfStart  [4097];

    // translation finite state machine base codes - ncbi4na
    enum EBaseCode {
        eBase_gap = 0,
        eBase_A,      /* A    */
        eBase_C,      /* C    */
        eBase_M,      /* AC   */
        eBase_G,      /* G    */
        eBase_R,      /* AG   */
        eBase_S,      /* CG   */
        eBase_V,      /* ACG  */
        eBase_T,      /* T    */
        eBase_W,      /* AT   */
        eBase_Y,      /* CT   */
        eBase_H,      /* ACT  */
        eBase_K,      /* GT   */
        eBase_D,      /* AGT  */
        eBase_B,      /* CGT  */
        eBase_N       /* ACGT */
    };

// initialize genetic code specific translation tables
    static void x_InitFsaTransl (std::string ncbieaa,
                                 std::string sncbieaa)
    {
        char        ch, aa, orf;
        bool        go_on;
        int         i, j, k, p, q, r, x, y, z, st, cd;
        static int  expansions [4] = {eBase_A, eBase_C, eBase_G, eBase_T};
        // T = 0, C = 1, A = 2, G = 3
        static int  codonIdx [9] = {0, 2, 1, 0, 3, 0, 0, 0, 0};

        // return if unable to find ncbieaa and sncbieaa strings
        if (ncbieaa == 0 || sncbieaa == 0) return;

        // also check length of ncbieaa and sncbieaa strings
        if (ncbieaa->size () != 64 || sncbieaa->size () != 64) return;

        // ambiguous codons map to unknown amino acid or not start
        for (i = 0; i <= 4096; i++) {
            m_AminoAcid [i] = 'X';
            m_OrfStart [i] = '-';
        }

        // lookup amino acid for each codon in genetic code table
        for (i = eBase_gap, st = 1; i <= eBase_N; i++) {
            for (j = eBase_gap; j <= eBase_N; j++) {
                for (k = eBase_gap; k <= eBase_N; k++, st++) {
                    aa = '\0';
                    orf = '\0';
                    go_on = true;

                    // expand ambiguous IJK nucleotide symbols into component bases XYZ
                    for (p = 0; p < 4 && go_on; p++) {
                        x = expansions [p];
                        if ((x & i) != 0) {
                            for (q = 0; q < 4 && go_on; q++) {
                                y = expansions [q];
                                if ((y & j) != 0) {
                                    for (r = 0; r < 4 && go_on; r++) {
                                        z = expansions [r];
                                        if ((z & k) != 0) {

                                            // calculate offset in genetic code string

                                            // the T = 0, C = 1, A = 2, G = 3 order is
                                            // necessary because the genetic code strings
                                            // are presented in TCAG order in printed tables
                                            // and in the genetic code strings
                                            cd = 16 * codonIdx [x] + 4 * codonIdx [y] + codonIdx [z];

                                            // lookup amino acid for codon XYZ
                                            ch = (*ncbieaa) [cd];
                                            if (aa == '\0') {
                                                aa = ch;
                                            } else if (aa != ch) {
                                                // allow Asx (Asp or Asn) and Glx (Glu or Gln)
                                                if ((aa == 'B' || aa == 'D' || aa == 'N') &&
                                                    (ch == 'D' || ch == 'N')) {
                                                    aa = 'B';
                                                } else if ((aa == 'Z' || aa == 'E' || aa == 'Q') &&
                                                           (ch == 'E' || ch == 'Q')) {
                                                    aa = 'Z';
                                                } else if ((aa == 'J' || aa == 'I' || aa == 'L') &&
                                                           (ch == 'I' || ch == 'L')) {
                                                    aa = 'J';
                                                } else {
                                                    aa = 'X';
                                                }
                                            }

                                            // lookup translation start flag
                                            ch = (*sncbieaa) [cd];
                                            if (orf == '\0') {
                                                orf = ch;
                                            } else if (orf != ch) {
                                                orf = 'X';
                                            }

                                            // drop out of loop as soon as answer is known
                                            if (aa == 'X' && orf == 'X') {
                                                go_on = false;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // assign amino acid and orf start
                    if (aa != '\0') {
                        m_AminoAcid [st] = aa;
                    }
                    if (orf != '\0') {
                        m_OrfStart [st] = orf;
                    }
                }
            }
        }
    }

};


#endif //MMSEQS_TRANSLATE_H
