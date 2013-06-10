#ifndef _SCORING_H
#define _SCORING_H 

#include <iostream>
#include <string>
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <malloc.h>
#include <xmmintrin.h>

#include "params.h"

// returns the array for the conversion of the amino acid number to aa one letter code
char* get_int2aa();

// returns the array for the conversion of the aa one letter code to the amino acid number (inversion of int2aa array)
int* get_aa2int(char* int2aa);

// reads the scoring matrix from the file
// (.mat files)
int** getScoringMatrix(std::string mFile, int* aa2int);

// reads the unbiased scoring matrix from the file containing the amino acid probabilities and returns the biased normalized score
// (.out files)
int** getBiasedScoringMatrix(std::string mFile, int* aa2int);

// get the smallest value in the scoring matrix that is used as bias in the scoring profile (s. Farrar paper section 2.2.2.)
unsigned char getMatrixMinValue(int** sc_matrix);

// generates the query profile from the query sequence and the scoring matrix (char values, range 0 <= score <= 255)
unsigned char* getQueryProfileByte(unsigned char* query, int qlen, int** sc_matrix, unsigned char bias);

// generates the query profile from the query sequence and the scoring matrix (short integer values, range 0 <= score <= 65535)
unsigned short* getQueryProfileWord(unsigned char* query, int qlen, int** sc_matrix);

// log2 calculation
inline float _log2 (float x) { return log10(x)/0.301029996; }

#endif
