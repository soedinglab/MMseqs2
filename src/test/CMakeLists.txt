include(MMseqsSetupTest)

set(TESTS
        #TestAdjustedKmerIterator.cpp
        TestAlignment.cpp
        TestAlignmentPerformance.cpp
        TestAlignmentTraceback.cpp
        TestAlp.cpp
        TestBacktraceTranslator.cpp
        TestCompositionBias.cpp
        TestCounting.cpp
        TestDBReader.cpp
        TestDBReaderIndexSerialization.cpp
        TestDiagonalScoring.cpp
        TestDiagonalScoringPerformance.cpp
        TestKmerGenerator.cpp
        TestKmerNucl.cpp
        TestKmerScore.cpp
        TestKwayMerge.cpp
        TestMultipleAlignment.cpp
        TestProfileAlignment.cpp
        TestPSSM.cpp
        TestPSSMPrune.cpp
        TestDBReaderZstd.cpp
        TestReduceMatrix.cpp
        TestScoreMatrixSerialization.cpp
        TestSequenceIndex.cpp
        TestTanTan.cpp
        TestTaxonomy.cpp
        TestTranslate.cpp
        TestTinyExpr.cpp
        TestTaxExpr.cpp
        TestUtil.cpp
        TestKsw2.cpp
        TestBestAlphabet.cpp
        )


FOREACH (TEST ${TESTS})
    mmseqs_setup_test(${TEST})
ENDFOREACH ()
