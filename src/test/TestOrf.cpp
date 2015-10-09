#include "gmock/gmock.h"

#include <Orf.h>
#include <vector>

class OrfTest : public testing::Test {
protected:

    const char* sequence = "CGAAGCGGGTGATGGCCGGCGCCGCGCCGGTTGGCGGCTGGCCATTCAAGGAGTGAGGAGATGGTCACTGGGCAGCGCGCCGGGGGGCGGCAGCAGCCCAAGGGTCGGGTCATTCCCGATTGGCCGCACCAGGCGCCCGCCACAGCCGGA";

    const char* reverseComplement = "TCCGGCTGTGGCGGGCGCCTGGTGCGGCCAATCGGGAATGACCCGACCCTTGGGCTGCTGCCGCCCCCCGGCGCGCTGCCCAGTGACCATCTCCTCACTCCTTGAATGGCCAGCCGCCAACCGGCGCGGCGCCGGCCATCACCCGCTTCG";

    size_t sequenceLength = 150;

    Orf* orf;

    virtual void SetUp() {
        orf = new Orf(sequence);
    }

    virtual void TearDown() {
        delete orf;
    }
};

TEST_F(OrfTest, Frame_1) {
    std::vector<std::string> computedOrfs;

    std::vector<Orf::SequenceLocation> results;
    FindForwardOrfs(sequence, sequenceLength, results, 1, 300, 0, Orf::FRAME_1, Orf::STRAND_PLUS);
    for(std::vector<Orf::SequenceLocation>::const_iterator it = results.begin(); it != results.end(); it++) {
        Orf::SequenceLocation loc = *it;
        auto buffer = orf->View(loc);
        computedOrfs.emplace_back(buffer.get());
    }

    std::vector<std::string> expectedOrfs;
    expectedOrfs.emplace_back("CGAAGCGGGTGA");
    expectedOrfs.emplace_back(
        "ATGGTCACTGGGCAGCGCGCCGGGGGGCGGCAGCAGCCCAAGGGTCGGGTCATTCCCGATTGGCCGCACCAGGCGCCCGCCACAGCCGGA"
    );

    EXPECT_THAT(computedOrfs, testing::ContainerEq(expectedOrfs));
}

TEST_F(OrfTest, Frame_2) {
    std::vector<std::string> computedOrfs;

    std::vector<Orf::SequenceLocation> results;
    FindForwardOrfs(sequence, sequenceLength, results, 1, 300, 0, Orf::FRAME_2, Orf::STRAND_PLUS);
    for(std::vector<Orf::SequenceLocation>::const_iterator it = results.begin(); it != results.end(); it++) {
        Orf::SequenceLocation loc = *it;
        auto buffer = orf->View(loc);\
        computedOrfs.emplace_back(buffer.get());
    }

    std::vector<std::string> expectedOrfs;
    expectedOrfs.emplace_back(
"GAAGCGGGTGATGGCCGGCGCCGCGCCGGTTGGCGGCTGGCCATTCAAGGAGTGAGGAGATGGTCACTGGGCAGCGCGCCGGGGGGCGGCAGCAGCCCAAGGGTCGGGTCATTCCCGATTGGCCGCACCAGGCGCCCGCCACAGCCGGA"
    );

    EXPECT_THAT(computedOrfs, testing::ContainerEq(expectedOrfs));
}

TEST_F(OrfTest, Frame_3) {
    std::vector<std::string> computedOrfs;

    std::vector<Orf::SequenceLocation> results;
    FindForwardOrfs(sequence, sequenceLength, results, 1, 300, 0, Orf::FRAME_3, Orf::STRAND_PLUS);
    for(std::vector<Orf::SequenceLocation>::const_iterator it = results.begin(); it != results.end(); it++) {
        Orf::SequenceLocation loc = *it;
        auto buffer = orf->View(loc);
        computedOrfs.emplace_back(buffer.get());
    }

    std::vector<std::string> expectedOrfs;
    expectedOrfs.emplace_back("ATGGCCGGCGCCGCGCCGGTTGGCGGCTGGCCATTCAAGGAGTGA");

    EXPECT_THAT(computedOrfs, testing::ContainerEq(expectedOrfs));
}

TEST_F(OrfTest, Frame_R_1) {
    std::vector<std::string> computedOrfs;

    std::vector<Orf::SequenceLocation> results;
    FindForwardOrfs(reverseComplement, sequenceLength, results, 1, 300, 0, Orf::FRAME_1, Orf::STRAND_MINUS);
    for(std::vector<Orf::SequenceLocation>::const_iterator it = results.begin(); it != results.end(); it++) {
        Orf::SequenceLocation loc = *it;
        auto buffer = orf->View(loc);
        computedOrfs.emplace_back(buffer.get());
    }

    std::vector<std::string> expectedOrfs;
    expectedOrfs.emplace_back(
"TCCGGCTGTGGCGGGCGCCTGGTGCGGCCAATCGGGAATGACCCGACCCTTGGGCTGCTGCCGCCCCCCGGCGCGCTGCCCAGTGACCATCTCCTCACTCCTTGA"
    );
    expectedOrfs.emplace_back("ATGGCCAGCCGCCAACCGGCGCGGCGCCGGCCATCACCCGCTTCG");

    EXPECT_THAT(computedOrfs, testing::ContainerEq(expectedOrfs));
}

TEST_F(OrfTest, Frame_R_2) {
    std::vector<std::string> computedOrfs;

    std::vector<Orf::SequenceLocation> results;
    FindForwardOrfs(reverseComplement, sequenceLength, results, 1, 300, 0, Orf::FRAME_3, Orf::STRAND_MINUS);
    for(std::vector<Orf::SequenceLocation>::const_iterator it = results.begin(); it != results.end(); it++) {
        Orf::SequenceLocation loc = *it;
        auto buffer = orf->View(loc);
        computedOrfs.emplace_back(buffer.get());
    }

    std::vector<std::string> expectedOrfs;
    expectedOrfs.emplace_back("CGGCTGTGGCGGGCGCCTGGTGCGGCCAATCGGGAATGA");

    EXPECT_THAT(computedOrfs, testing::ContainerEq(expectedOrfs));
}

TEST_F(OrfTest, Frame_R_3) {
    std::vector<std::string> computedOrfs;

    std::vector<Orf::SequenceLocation> results;
    FindForwardOrfs(reverseComplement, sequenceLength, results, 1, 300, 0, Orf::FRAME_2, Orf::STRAND_MINUS);
    for(std::vector<Orf::SequenceLocation>::const_iterator it = results.begin(); it != results.end(); it++) {
        Orf::SequenceLocation loc = *it;
        auto buffer = orf->View(loc);
        computedOrfs.emplace_back(buffer.get());
    }

    std::vector<std::string> expectedOrfs;
    expectedOrfs.emplace_back(
"ATGACCCGACCCTTGGGCTGCTGCCGCCCCCCGGCGCGCTGCCCAGTGACCATCTCCTCACTCCTTGAATGGCCAGCCGCCAACCGGCGCGGCGCCGGCCATCACCCGCTTCG"
    );

    EXPECT_THAT(computedOrfs, testing::ContainerEq(expectedOrfs));
}

TEST_F(OrfTest, Orf_All) {
    std::vector<Orf::SequenceLocation> computedOrfs;

    orf->FindOrfs(computedOrfs);


    EXPECT_EQ(computedOrfs.size(), 8);
}