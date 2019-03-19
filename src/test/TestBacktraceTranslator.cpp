#include "BacktraceTranslator.h"
#include "Parameters.h"

const char* binary_name = "test_backtracetranslator";

int main(int, const char**) {
    // s1 5 ATT-GCA 11
    // s2 3 ATTTG-- 8
    // s3 6 --T-G-A 9
    // AB, BC

    // ATTTG-- MMMIM
    // ATT-GCA

    // ATTGCA MMMIM
    // --TG-A MMIM
    // s2 -> s1, s1-> s3 => infer s2 -> s1 -> s3
//    Matcher::result_t resultAB(2, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                           // qs qe      ts te
//                               3, 8, 10,  5, 9, 15, "MMMIM");
//    Matcher::result_t resultBC(3, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                               7, 11, 15, 6, 9, 20, "MMIM");


    // ATT-G-- MMMIM
    // ATTTGCA

    // ATTTGCA MMMIM
    // --T-GCA MIMMM
    // s2 -> s1, s1-> s3 => infer s2 -> s1 -> s3
//    Matcher::result_t resultAB(2, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                           // qs qe      ts te
//                               0, 3, 10,  0, 4, 15, "MMMDM");
//    Matcher::result_t resultBC(3, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                               2, 6, 15, 0, 3, 20, "MIMMM");



//    // ATTTTG-- MMMIM
//    // ATT-TGCA
//
//    // ATT-TGCA MMMIM
//    // --TTTGC- MIMMM
//    // s2 -> s1, s1-> s3 => infer s2 -> s1 -> s3
//    Matcher::result_t resultAB(2, 8, 0.6, 0.8, 0.8, 0.001, 6,
////                           // qs qe      ts te
//                               0, 5, 10,  0, 3, 15, "MMMIMM");
//    Matcher::result_t resultBC(3, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                               2, 5, 15, 0, 4, 20, "MDMMM");


    // ATT-TG-- MMMIM
    // ATTTTGCA

    // ATT-TGCA MMMIM
    // --TTTGC- MIMMM
    // s2 -> s1, s1-> s3 => infer s2 -> s1 -> s3
    Matcher::result_t resultAB(2, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                           // qs qe      ts te
                               0, 4, 10,  0, 5, 15, "MMMDMM");
    Matcher::result_t resultBC(3, 8, 0.6, 0.8, 0.8, 0.001, 6,
                               2, 5, 15, 0, 4, 20, "MDMMM");

    // ATT-G-- MMMIM
    // ATTTGCA

    // ATTTGCA MMMIM
    // --T-GCA MIMMM
    // s2 -> s1, s1-> s3 => infer s2 -> s1 -> s3
//    Matcher::result_t resultAB(2, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                           // qs qe      ts te
//                               0, 3, 10,  0, 4, 15, "MMMDM");
//    Matcher::result_t resultBC(3, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                               2, 6, 15, 0, 3, 20, "MIMMM");

//
//    Matcher::result_t resultAC;
//    resultAC.backtrace.reserve(2000);
//    BacktraceTranslator translator;
//    translator.translateResult(resultAB, resultBC, resultAC);
//


    // s1 5 AT--GCA 11
    // s2 3 ATTTG-- 8
    // s3 6 --T-G-A 9
    // AB, BC

    // AT---G MMMDDM
    // ATTTAG

    // ATTTAG
    // -TTT-G MMIM

    // ATTG
    // --TG
    // s2 -> s1, s1-> s3 => infer s2 -> s1 -> s3
//    Matcher::result_t resultAB(2, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                           // qs qe      ts te
//                               0, 3, 10,  0, 5, 15, "MMDDDM");
//    Matcher::result_t resultBC(3, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                               1, 5, 15, 0, 4, 20, "MMMIM");


    // s1 5 AT--GCA 11
    // s2 3 ATTTG-- 8
    // s3 6 --T-G-A 9
    // AB, BC

    // AT---G MMMDDM
    // ATTTAG

    // ATTTAG
    // --TT-- MMIM

    // ATTG
    // --TG
    // s2 -> s1, s1-> s3 => infer s2 -> s1 -> s3
//    Matcher::result_t resultAB(2, 8, 0.6, 0.8, 0.8, 0.001, 6,
//             qs qe      ts te
//                               0, 3, 10,  0, 5, 15, "MMDDDM");
//    Matcher::result_t resultBC(3, 8, 0.6, 0.8, 0.8, 0.001, 6,
//                               2, 5, 15, 0, 2, 20, "MM");
//    MM
//    IM
//    DM
//    MD
//    DD
//    MI
//    II
//    DI


    Matcher::result_t resultAC;
    resultAC.backtrace.reserve(2000);
    BacktraceTranslator translator;
    translator.translateResult(resultAB, resultBC, resultAC);

    //  0| 2|
    //   ELLO 4
    // SHELL
    // 1| 4| 5

    char buffer[2048];
    Matcher::resultToBuffer(buffer, resultAC, true, true);
    Debug(Debug::INFO) << buffer;

    return EXIT_SUCCESS;
}
