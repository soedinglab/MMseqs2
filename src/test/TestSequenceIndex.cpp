#include <iostream>
#include <list>
#include <algorithm>
#include <cmath>

#include "SequenceLookup.h"
#include "SubstitutionMatrix.h"
#include "Clustering.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Parameters.h"

const char* binary_name = "test_sequenceindex";
DEFAULT_PARAMETER_SINGLETON_INIT

int main (int, const char**) {
    size_t kmer_size = 6;
    Parameters &par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    std::string S1 = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF";
    const char *S1char = S1.c_str();
    std::cout << S1char << "\n\n";
    Sequence s1(10000, 0, &subMat, kmer_size, true, false);
    s1.mapSequence(0, 0, S1char, strlen(S1char));
    std::string S2 = "PQFSLWKRPVVTAYIEGQPVEVLLDTGADDSIVAGIELGNNIVGGIGGFINTLEYKNVEIEVLNKKVRATIMTGDTPINIFGRNILTALGMSLNL";
    const char *S2char = S2.c_str();
    std::cout << S2char << "\n\n";
    Sequence s2(10000,  0, &subMat, kmer_size, true, false);
    s2.mapSequence(1, 1, S2char, strlen(S2char));
    std::string S3 = "PQFHLWKRPVVTAGQPVEVLLDTGADDSIVTGIELGPHYTPKIVGGIGGFINTKEYKNVEVEVLGKRIKGTIMTGDTPINIFGRNLLTALGMSLNF";
    const char *S3char = S3.c_str();
    std::cout << S3char << "\n\n";
    Sequence s3(10000,  0, &subMat, kmer_size, true, false);
    s3.mapSequence(2, 2, S3char, strlen(S3char));
    std::string S4 = "LAMTMEHKDRPLVRVILTNTGSHPVKQRSVYITALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL";
    const char *S4char = S4.c_str();
    std::cout << S4char << "\n\n";
    Sequence s4(10000,  0, &subMat, kmer_size, true, false);
    s4.mapSequence(3, 3, S4char, strlen(S4char));

    SequenceLookup lookup(4, s1.L + s2.L + s3.L + s4.L);
    lookup.addSequence(&s1);
    lookup.addSequence(&s2);
    lookup.addSequence(&s3);
    lookup.addSequence(&s4);

    std::pair<const unsigned char *, const unsigned int> s1res = lookup.getSequence(0);
    std::cout << (const char *) s1res.first << std::endl;
    if (s1res.second != S1.length())
        std::cout << "Diff length" << std::endl;
    for (size_t i = 0; i < s1res.second; i++) {
        if (subMat.num2aa[s1res.first[i]] != S1char[i]) {
            std::cout << "Wrong data" << std::endl;
        }
    }
    std::pair<const unsigned char *, const unsigned int> s2res = lookup.getSequence(1);
    if (s2res.second != S2.length())
        std::cout << "Diff length" << std::endl;
    for (size_t i = 0; i < s2res.second; i++) {
        if (subMat.num2aa[s2res.first[i]] != S2char[i]) {
            std::cout << "Wrong data" << std::endl;
        }
    }

    std::pair<const unsigned char *, const unsigned int> s3res = lookup.getSequence(2);
    if (s3res.second != S3.length())
        std::cout << "Diff length" << std::endl;

    for (size_t i = 0; i < s3res.second; i++) {
        if (subMat.num2aa[s3res.first[i]] != S3char[i]) {
            std::cout << "Wrong data" << std::endl;
        }
    }

    std::pair<const unsigned char *, const unsigned int> s4res = lookup.getSequence(3);
    for (size_t i = 0; i < s4res.second; i++) {
        if (subMat.num2aa[s4res.first[i]] != S4char[i]) {
            std::cout << "Wrong data" << std::endl;
        }
    }
}
