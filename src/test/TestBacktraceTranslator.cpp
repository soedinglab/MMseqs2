#include "BacktraceTranslator.h"
#include "Parameters.h"

const char* binary_name = "test_backtracetranslator";

int main(int, const char**) {
//    const char* dataAB = "76\t2285\t1.000\t0.000E+00\t0\t1123\t1124\t0\t1123\t1124\t1124M\n";
//    const char* dataBC = "43\t822\t0.524\t1.911E-259\t635\t1122\t1124\t7\t507\t575\t31M5D55M5D43M7D153M1I9M4I26M1I146M1D14M1D5M\n";
    //const char* dataBC = "3\t1184\t1.000\t0.000E+00\t0\t574\t575\t0\t574\t575\t575M\n";


//    const char* dataAB = "1\t0\t1.000\t0.000E+00\t0\t5\t6\t0\t5\t6\t5M\n";
//    const char* dataBC = "2\t0\t1.000\t0.000E+00\t0\t5\t6\t0\t5\t6\t5M\n";
    // 1| 4|
    // HELLO 5
    // JELLO
    // 1| 4| 5
//    const char* dataAB = "1\t0\t1.000\t0.000E+00\t1\t4\t5\t1\t4\t5\t4M\n";
    //  1|3|
    //  JELLO 5
    // SHELL
    //  2|4|  5
//     const char* dataBC = "2\t0\t1.000\t0.000E+00\t1\t4\t5\t2\t4\t5\t3M\n";

    //      2|4|
    //     JELLO 5
    // ARMADILLO 9
    //      6|8|
//    const char* dataBC = "2\t0\t1.000\t0.000E+00\t2\t4\t5\t6\t8\t9\t3M\n";


    // 0| 3|
    //  ELLO 4
    // HELLO
    // 1| 4| 5
    const char* dataAB = "1\t0\t1.000\t0.000E+00\t0\t3\t4\t1\t4\t5\t4M\n";
    // 0| 3|
    //  HELLO 5
    // SHELL
    // 1| 4| 5
    const char* dataBC = "2\t0\t1.000\t0.000E+00\t0\t3\t5\t1\t4\t5\t4M\n";
    Debug(Debug::INFO) << dataAB;
    Debug(Debug::INFO) << dataBC;

    Matcher::result_t resultAB = Matcher::parseAlignmentRecord(dataAB, false);
    Matcher::result_t resultBC = Matcher::parseAlignmentRecord(dataBC, false);

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
