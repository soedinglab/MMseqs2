///////////////////////////////////////////////////////////////
// Constants for the SSE2 Smith Waterman alignment calculation
///////////////////////////////////////////////////////////////

// size of amino acid alphabet
const int AMINOACID_DIM = 20;
// bias in the scoring matrix, makes alignments shorter (only used for the calculation of the biased matrix)
const float SCORING_BIAS = -0.2;
// gap open penalty
const unsigned char GAP_OPEN = 11;
// gap extension penalty
const unsigned char GAP_EXTEND = 1;
// Smith Waterman score threshold
const float SW_SCORE_THR = 0;
// Bit factor in the scoring matrix
const double BIT_FACTOR = 2.0;
