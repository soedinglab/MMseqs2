include(MMseqsSetupTest)

set(TESTS
        TestAlignment.cpp
        TestAlignmentPerformance.cpp
        TestAlignmentTraceback.cpp
        TestAlp.cpp
        TestCompositionBias.cpp
        TestCounting.cpp
        TestDBReader.cpp
        TestDBReaderIndexSerialization.cpp
        TestDiagonalScoring.cpp
        TestDiagonalScoringPerformance.cpp
        TestIndexTable.cpp
        TestKmerGenerator.cpp
        TestKmerScore.cpp
        TestMultipleAlignment.cpp
        TestProfileAlignment.cpp
        TestPSSM.cpp
        TestPSSMPrune.cpp
        TestReduceMatrix.cpp
        TestScoreMatrixSerialization.cpp
        TestSequenceIndex.cpp
        TestTanTan.cpp
        TestTranslate.cpp
        TestUtil.cpp
        )

FOREACH (TEST ${TESTS})
    mmseqs_setup_test(${TEST})
ENDFOREACH ()
