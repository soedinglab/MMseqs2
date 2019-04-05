/* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*/

#ifndef MMSEQS_TRANSLATE_H
#define MMSEQS_TRANSLATE_H

#include <string>
#include "Debug.h"
#include "Util.h"
#include <set>
#include <cmath>

// standard genetic code
//
// ncbieaa  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
// sncbieaa "---M---------------M---------------M----------------------------"
//
// -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
// -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
// -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

/*
                         Second Position
First      T             C             A             G      Third
-----------------------------------------------------------------
  T   TTT Phe [F]   TCT Ser [S]   TAT Tyr [Y]   TGT Cys [C]   T
      TTC Phe [F]   TCC Ser [S]   TAC Tyr [Y]   TGC Cys [C]   C
      TTA Leu [L]   TCA Ser [S]   TAA Ter [*]   TGA Ter [*]   A
      TTG Leu [L]   TCG Ser [S]   TAG Ter [*]   TGG Trp [W]   G
-----------------------------------------------------------------
  C   CTT Leu [L]   CCT Pro [P]   CAT His [H]   CGT Arg [R]   T
      CTC Leu [L]   CCC Pro [P]   CAC His [H]   CGC Arg [R]   C
      CTA Leu [L]   CCA Pro [P]   CAA Gln [Q]   CGA Arg [R]   A
      CTG Leu [L]   CCG Pro [P]   CAG Gln [Q]   CGG Arg [R]   G
-----------------------------------------------------------------
  A   ATT Ile [I]   ACT Thr [T]   AAT Asn [N]   AGT Ser [S]   T
      ATC Ile [I]   ACC Thr [T]   AAC Asn [N]   AGC Ser [S]   C
      ATA Ile [I]   ACA Thr [T]   AAA Lys [K]   AGA Arg [R]   A
      ATG Met [M]   ACG Thr [T]   AAG Lys [K]   AGG Arg [R]   G
-----------------------------------------------------------------
  G   GTT Val [V]   GCT Ala [A]   GAT Asp [D]   GGT Gly [G]   T
      GTC Val [V]   GCC Ala [A]   GAC Asp [D]   GGC Gly [G]   C
      GTA Val [V]   GCA Ala [A]   GAA Glu [E]   GGA Gly [G]   A
      GTG Val [V]   GCG Ala [A]   GAG Glu [E]   GGG Gly [G]   G
-----------------------------------------------------------------
*/

// local copy of gc.prt genetic code table ASN.1
//const char * CGen_code_table_imp::sm_GenCodeTblMemStr [] =
//        {
//                "Genetic-code-table ::= {\n",
// CANONICAL
//                "{ name \"Standard\" , name \"SGC0\" , id 1 ,\n",
//                "ncbieaa  \"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",\n",
//                "sncbieaa \"---M------**--*----M---------------M----------------------------\" } ,\n",


class TranslateNucl {
public:
    enum GenCode {
        CANONICAL = 1,
        VERT_MITOCHONDRIAL,
        YEAST_MITOCHONDRIAL,
        MOLD_MITOCHONDRIAL,
        INVERT_MITOCHONDRIAL,
        CILIATE,
        FLATWORM_MITOCHONDRIAL = 9,
        EUPLOTID,
        PROKARYOTE,
        ALT_YEAST,
        ASCIDIAN_MITOCHONDRIAL,
        ALT_FLATWORM_MITOCHONDRIAL,
        BLEPHARISMA,
        CHLOROPHYCEAN_MITOCHONDRIAL,
        TREMATODE_MITOCHONDRIAL = 21,
        SCENEDESMUS_MITOCHONDRIAL,
        THRAUSTOCHYTRIUM_MITOCHONDRIAL,
        PTEROBRANCHIA_MITOCHONDRIAL,
        GRACILIBACTERIA,
        PACHYSOLEN,
        KARYORELICT,
        CONDYLOSTOMA,
        MESODINIUM,
        PERTRICH,
        BLASTOCRITHIDIA
    };

    TranslateNucl(GenCode code){
        std::string ncbieaa = "";
        std::string sncbieaa = "";
        switch(code){
            case CANONICAL:
//                "{ name \"Standard\" , name \"SGC0\" , id 1 ,\n",
                ncbieaa ="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="---M------**--*----M---------------M----------------------------";
                break;
            case VERT_MITOCHONDRIAL:
//                "{ name \"Vertebrate Mitochondrial\" , name \"SGC1\" , id 2 ,\n",
                ncbieaa ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
                sncbieaa="----------**--------------------MMMM----------**---M------------";
                break;
            case YEAST_MITOCHONDRIAL:
//                "{ name \"Yeast Mitochondrial\" , name \"SGC2\" , id 3 ,\n",
                ncbieaa ="FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="----------**----------------------MM----------------------------";
                break;
            case MOLD_MITOCHONDRIAL:
//                "{ name \"Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate\n",
//                "Mitochondrial; Mycoplasma; Spiroplasma\" ,\n",
//                "name \"SGC3\" , id 4 ,\n",
                ncbieaa ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="--MM------**-------M------------MMMM---------------M------------";
                break;
            case INVERT_MITOCHONDRIAL:
//                "{ name \"Invertebrate Mitochondrial\" , name \"SGC4\" , id 5 ,\n",
                ncbieaa ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
                sncbieaa="---M------**--------------------MMMM---------------M------------";
                break;
            case CILIATE:
//                "{ name \"Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear\" ,\n",
//                "name \"SGC5\" , id 6 ,\n",
                ncbieaa ="FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="--------------*--------------------M----------------------------";
                break;
            case FLATWORM_MITOCHONDRIAL:
//                "{ name \"Echinoderm Mitochondrial; Flatworm Mitochondrial\" , name \"SGC8\" , id 9 ,\n",
                ncbieaa ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
                sncbieaa="----------**-----------------------M---------------M------------";
                break;
            case EUPLOTID:
//                "{ name \"Euplotid Nuclear\" , name \"SGC9\" , id 10 ,\n",
                ncbieaa ="FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="----------**-----------------------M----------------------------";
                break;
            case PROKARYOTE:
//                "{ name \"Bacterial, Archaeal and Plant Plastid\" , id 11 ,\n",
                ncbieaa ="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="---M------**--*----M------------MMMM---------------M------------";
                break;
            case ALT_YEAST:
//                "{ name \"Alternative Yeast Nuclear\" , id 12 ,\n",
                ncbieaa ="FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="----------**--*----M---------------M----------------------------";
                break;
            case ASCIDIAN_MITOCHONDRIAL:
//                "{ name \"Ascidian Mitochondrial\" , id 13 ,\n",
                ncbieaa ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
                sncbieaa="---M------**----------------------MM---------------M------------";
                break;
            case ALT_FLATWORM_MITOCHONDRIAL:
//                "{ name \"Alternative Flatworm Mitochondrial\" , id 14 ,\n",
                ncbieaa ="FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
                sncbieaa="-----------*-----------------------M----------------------------";
                break;
            case BLEPHARISMA:
//                "{ name \"Blepharisma Macronuclear\" , id 15 ,\n",
                ncbieaa ="FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="----------*---*--------------------M----------------------------";
                break;
            case CHLOROPHYCEAN_MITOCHONDRIAL:
//                "{ name \"Chlorophycean Mitochondrial\" , id 16 ,\n",
                ncbieaa ="FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="----------*---*--------------------M----------------------------";
                break;
            case TREMATODE_MITOCHONDRIAL:
//                "{ name \"Trematode Mitochondrial\" , id 21 ,\n",
                ncbieaa ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
                sncbieaa="----------**-----------------------M---------------M------------";
                break;
            case SCENEDESMUS_MITOCHONDRIAL:
//                "{ name \"Scenedesmus obliquus Mitochondrial\" , id 22 ,\n",
                ncbieaa ="FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="------*---*---*--------------------M----------------------------";
                break;
            case THRAUSTOCHYTRIUM_MITOCHONDRIAL:
//                "{ name \"Thraustochytrium Mitochondrial\" , id 23 ,\n",
                ncbieaa ="FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="--*-------**--*-----------------M--M---------------M------------";
                break;
            case PTEROBRANCHIA_MITOCHONDRIAL:
//                "{ name \"Pterobranchia Mitochondrial\" , id 24 ,\n",
                ncbieaa ="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG";
                sncbieaa="---M------**-------M---------------M---------------M------------";
                break;
            case GRACILIBACTERIA:
//                "{ name \"Candidate Division SR1 and Gracilibacteria\" , id 25 ,\n",
                ncbieaa ="FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="---M------**-----------------------M---------------M------------";
                break;
            case PACHYSOLEN:
//                "{ name \"Pachysolen tannophilus Nuclear\" , id 26 ,\n",
                ncbieaa ="FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="----------**--*----M---------------M----------------------------";
                break;
            case KARYORELICT:
//                "{ name \"Karyorelict Nuclear\" , id 27 ,\n",
                ncbieaa ="FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="--------------*--------------------M----------------------------";
                break;
            case CONDYLOSTOMA:
//                "{ name \"Condylostoma Nuclear\" , id 28 ,\n",
                ncbieaa ="FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="----------**--*--------------------M----------------------------";
                break;
            case MESODINIUM:
//                "{ name \"Mesodinium Nuclear\" , id 29 ,\n",
                ncbieaa ="FFLLSSSSYYYYCC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="--------------*--------------------M----------------------------";
                break;
            case PERTRICH:
//                "{ name \"Peritrich Nuclear\" , id 30 ,\n",
                ncbieaa ="FFLLSSSSYYEECC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="--------------*--------------------M----------------------------";
                break;
            case BLASTOCRITHIDIA:
//                "{ name \"Blastocrithidia Nuclear\" , id 31 ,\n",
                ncbieaa ="FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
                sncbieaa="----------**-----------------------M----------------------------";
                break;
            default:
                Debug(Debug::ERROR) << "Invalid translation table selected!\n";
                EXIT(EXIT_FAILURE);
                break;
        }
        // init table
        initTranslationTable(&ncbieaa ,&sncbieaa);
        initConversionTable();
    };
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

    // set because we want unique keys. several states map to the same start/stop codon
    std::set<int> stopCodons;
    std::set<int> startCodons;

    std::vector<std::string> getStartCodons() {
        return getCodons(startCodons);
    }

    std::vector<std::string> getStopCodons() {
        return getCodons(stopCodons);
    }
    
    std::vector<std::string> getCodons(const std::set<int> &codonsSet) {        
        std::vector<std::string> codonsVec;
        for (std::set<int>::const_iterator it=codonsSet.begin(); it!=codonsSet.end(); ++it) {
            int currCode = *it;
            std::string currStr;
            for (size_t nucInd = 0; nucInd < 3; ++nucInd) {
                int currPower = (int) std::pow(4,(2 - nucInd));
                int currQ = currCode / currPower;

                // T = 0, C = 1, A = 2, G = 3
                char currNuc = 'T';
                if (currQ == 1) {
                    currNuc = 'C';
                }
                else if (currQ == 2) {
                    currNuc = 'A';
                }
                else if (currQ == 3) {
                    currNuc = 'G';
                }
                currStr.push_back(currNuc);

                int currR = currCode % currPower;
                currCode = currR;
            }
            codonsVec.push_back(currStr);
        }
        return (codonsVec);
    }

    // static instances of single copy translation tables common to all genetic codes
    int sm_NextState  [4097];
    int sm_RvCmpState [4097];
    int sm_BaseToIdx  [256];


    int getCodonState (int state, unsigned char ch) const {
        if (state < 0 || state > 4096) return 0;
        return (sm_NextState [state] + sm_BaseToIdx [(int) ch]);
    }

    char getCodonResidue (int state) const {
        if (state < 0 || state > 4096) return 0;
        return (m_AminoAcid[state]);
    }

    // initialize base conversion, next state, and reverse complement state tables
    void initConversionTable (void)
    {
        char         ch;
        int          i, j, k, p, q, r, nx, st;
        static char  charToBase [17] = "-ACMGRSVTWYHKDBN";
        static char  baseToComp [17] = "-TGKCYSBAWRDMHVN";

        // illegal characters map to 0
        for (i = 0; i < 256; i++) {
            sm_BaseToIdx [i] = 0;
        }

        // map iupacna alphabet to EBaseCode
        for (i = eBase_gap; i <= eBase_N; i++) {
            ch = charToBase [i];
            sm_BaseToIdx [(int) ch] = i;
            ch = (unsigned char)tolower (ch);
            sm_BaseToIdx [(int) ch] = i;
        }
        sm_BaseToIdx [(int) 'U'] = eBase_T;
        sm_BaseToIdx [(int) 'u'] = eBase_T;
        sm_BaseToIdx [(int) 'X'] = eBase_N;
        sm_BaseToIdx [(int) 'x'] = eBase_N;

        // also map ncbi4na alphabet to EBaseCode
        for (i = eBase_gap; i <= eBase_N; i++) {
            sm_BaseToIdx [(int) i] = i;
        }

        // treat state 0 as already having seen NN,
        // avoiding single and double letter states
        sm_NextState [0] = 4081;
        sm_RvCmpState [0] = 4096;

        // states 1 through 4096 are triple letter states (---, --A, ..., NNT, NNN)
        for (i = eBase_gap, st = 1; i <= eBase_N; i++) {
            for (j = eBase_gap, nx = 1; j <= eBase_N; j++) {
                for (k = eBase_gap; k <= eBase_N; k++, st++, nx += 16) {
                    sm_NextState [st] = nx;
                    p = sm_BaseToIdx [(int)  baseToComp [k]];
                    q = sm_BaseToIdx [(int)  baseToComp [j]];
                    r = sm_BaseToIdx [(int)  baseToComp [i]];
                    sm_RvCmpState [st] = 256 * p + 16 * q + r + 1;
                }
            }
        }
    }

// initialize genetic code specific translation tables
    void initTranslationTable (std::string * ncbieaa,
                               std::string * sncbieaa)
    {
        char        ch, aa, orf;
        bool        go_on;
        int         i, j, k, p, q, r, x, y, z, st, cd;
        static int  expansions [4] = {eBase_A, eBase_C, eBase_G, eBase_T};
        // T = 0, C = 1, A = 2, G = 3
        static int  codonIdx [9] = {0, 2, 1, 0, 3, 0, 0, 0, 0};
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
                                            ch = ncbieaa->at(cd);
                                            
                                            if (aa == '\0') {
                                                aa = ch;
                                                // here is a stop codon
                                                if (aa == '*') {
                                                    stopCodons.insert(cd);
                                                }
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
                                            ch = sncbieaa->at(cd);
                                            if (orf == '\0') {
                                                orf = ch;
                                            } else if (orf != ch) {
                                                orf = 'X';
                                            }
                                            // here is a start codon
                                            if (ch == 'M') {
                                                startCodons.insert(cd);
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
    };

    void translate(char *aa, const char *nucl, int L) const {
        int state = 0;
        for (int i = 0;  i < L;  i += 3) {
            // loop through one codon at a time
            bool isLowerCase = false;
            for (int k = 0;  k < 3;  ++k) {
                isLowerCase |= islower(nucl[i+k]);
                state = getCodonState(state, nucl[i+k]);
            }
//            std::cout << state  << " ";
            int pos = i/3;
            char residue = getCodonResidue(state);
            aa[pos] = (isLowerCase) ? tolower(residue) : residue;
        }
//        std::cout << std::endl;
    }

    char translateSingleCodon(const char *nucl) const {
        int state = 0;
        for (int k = 0;  k < 3;  ++k) {
            state = getCodonState(state, nucl[k]);
        }
        return getCodonResidue(state);
    }
};

#endif //MMSEQS_TRANSLATE_H
