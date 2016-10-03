//
// Created by lars on 26.05.15.
//

#ifndef MMSEQS_DISTANCECALCULATOR_H
#define MMSEQS_DISTANCECALCULATOR_H


#include "Debug.h"
#include <sstream>

class DistanceCalculator {
public:
    /*
     * Adapted from levenshtein.js (https://gist.github.com/andrei-m/982927)
     * Changed to hold only two rows of the dynamic programing matrix
     * Copyright (c) 2011 Andrei Mackenzie
     * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
     * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
     * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
    static size_t levenshteinDistance(const std::string &s1, const std::string &s2) {
        size_t m = s1.length();
        size_t n = s2.length();

        if (m == 0)
            return n;
        if (n == 0)
            return m;

        size_t *matrix[2];
        // initialize first row and first column with increasing gap costs
        for (size_t i = 0; i < 2; ++i) {
            matrix[i] = new size_t[m + 1];
            matrix[i][0] = i;
        }

        for (size_t i = 0; i <= m; ++i) {
            matrix[0][i] = i;
        }

        // fill in the rest of the matrix
        for (size_t i = 1; i <= n; i++) {
            // previous row can be overwritten now
            size_t *curr = matrix[i & 1];
            size_t *prev = matrix[!(i & 1)];
            curr[0] = i;
            for (size_t j = 1; j <= m; j++) {
                if (s2[i - 1] == s1[j - 1]) {
                    curr[j] = prev[j - 1]; // match
                } else {
                    curr[j] = std::min(prev[j - 1] + 1, // substitution
                              std::min(curr[j - 1] + 1, // insertion
                                       prev[j] + 1));   // deletion
                }
            }
        }

        // the last element of the matrix is the distance
        size_t distance = matrix[n & 1][m];

        for (size_t i = 0; i < 2; ++i) {
            delete matrix[i];
        }

        return distance;
    }

    static double keywordDistance(const std::string &s1, const std::string &s2) {
        std::stringstream seq1(s1);
        std::string segment;
        std::vector<std::string> kwlist1;
        while (std::getline(seq1, segment, ';')) {
            if (strcmp(segment.c_str(), "\n") != 0)
                kwlist1.push_back(segment);
        }
        std::vector<std::string> kwlist2;
        std::stringstream seq2(s2);
        while (std::getline(seq2, segment, ';')) {
            if (strcmp(segment.c_str(), "\n") != 0)
                kwlist2.push_back(segment);
        }
        double score = 0.0;
        for (std::string s1 : kwlist1) {
            for (std::string s2 : kwlist2) {
                if (strcmp(s1.c_str(), s2.c_str()) == 0) {
                    score++;
                    break;
                }
            }
        }
        score = ((score / std::max(kwlist1.size(), kwlist2.size())) +
                 (score / std::min(kwlist1.size(), kwlist2.size()))) / 2;

        return score;
    }

};

#endif //MMSEQS_DISTANCECALCULATOR_H
