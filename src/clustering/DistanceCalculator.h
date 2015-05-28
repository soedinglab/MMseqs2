//
// Created by lars on 26.05.15.
//

#ifndef MMSEQS_DISTANCECALCULATOR_H
#define MMSEQS_DISTANCECALCULATOR_H


#include <Debug.h>
#include <sstream>

class DistanceCalculator {
public:

    // Compute Levenshtein Distance
// Martin Ettl, 2012-10-05
    static size_t uiLevenshteinDistance(const std::string &s1, const std::string &s2)
    {
        const size_t m(s1.size());
        const size_t n(s2.size());

        if( m==0 ) return n;
        if( n==0 ) return m;

        size_t *costs = new size_t[n + 1];

        for( size_t k=0; k<=n; k++ ) costs[k] = k;

        size_t i = 0;
        for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
        {
            costs[0] = i+1;
            size_t corner = i;

            size_t j = 0;
            for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
            {
                size_t upper = costs[j+1];
                if( *it1 == *it2 )
                {
                    costs[j+1] = corner;
                }
                else
                {
                    size_t t(upper<corner?upper:corner);
                    costs[j+1] = (costs[j]<t?costs[j]:t)+1;
                }

                corner = upper;
            }
        }

        size_t result = costs[n];
        delete [] costs;

        return result;
    }
    static double keywordDistance(const std::string &s1, const std::string &s2){
        std::stringstream seq1(s1);
        std::string segment;
        std::vector<std::string> kwlist1;
        while(std::getline(seq1, segment, ';'))
        {
            if(strcmp(segment.c_str(),"\n") != 0)
            kwlist1.push_back(segment);
        }
        std::vector<std::string> kwlist2;
        std::stringstream seq2(s2);
        while(std::getline(seq2, segment, ';'))
        {
            if(strcmp(segment.c_str(),"\n")!=0)
            kwlist2.push_back(segment);
        }
    double score=0.0;
        for(std::string s1 : kwlist1){
            for(std::string s2 : kwlist2){
                if(strcmp(s1.c_str(),s2.c_str())==0){
                    score++;
                    break;
                }
            }
        }
        score=((score/std::max(kwlist1.size(),kwlist2.size()))+(score/std::min(kwlist1.size(),kwlist2.size())))/2;

       // Debug(Debug::INFO) << s1 <<"compared to "<<s2 <<"score:"<<score <<"\n";
        return score;
    }

};

#endif //MMSEQS_DISTANCECALCULATOR_H
