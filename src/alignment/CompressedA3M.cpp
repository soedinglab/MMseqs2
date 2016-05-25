//
// Created by Milot on 24/05/16.
//

#include "CompressedA3M.h"
#include "DBConcat.h"
#include <sstream>

std::string CompressedA3M::fromAlignmentResult(const std::vector<Matcher::result_t>& alignment, DBConcat* referenceDBr) {
    std::ostringstream out;

    int totalNumberOfBlocks = 0;
    char zero = 0;
    // then write the alignment
    for (size_t i = 0; i < alignment.size(); i++) {
        unsigned int nbOfBlocks = 0;
        std::stringstream blocksDescription, firstBlock;

        Matcher::result_t aln = alignment.at(i);

        // detect the blocks
        for (size_t btIndex = 0; btIndex < (aln.backtrace).size();) {
            int indelLen = 0;
            int matchLen = 0;
            char inOrDel = 0;
            // seek the next insertion or deletion
            while (btIndex < (aln.backtrace).size() && aln.backtrace.at(btIndex) == 'M' &&
                   matchLen < 255) {
                btIndex++;
                matchLen++;
            }

            if (btIndex < (aln.backtrace).size() &&
                aln.backtrace.at(btIndex) != 'M') // end of block because an I or D was found
                inOrDel = aln.backtrace.at(btIndex); // store whether it is I or D

            // seek the next match
            while (btIndex < (aln.backtrace).size() && aln.backtrace.at(btIndex) == inOrDel &&
                   indelLen < 127) {
                btIndex++;
                indelLen++;
            }
            // deletion must be count negatively
            // and deletion states are insertion states from the sight of the aligned sequence
            if (indelLen && inOrDel == 'I')
                indelLen *= -1;

            blocksDescription.write(reinterpret_cast<const char *>(&matchLen), sizeof(char));
            blocksDescription.write(reinterpret_cast<const char *>(&indelLen), sizeof(char));
            nbOfBlocks++;
        }

        // add the deletion at the begining of the sequence (local alignment)
        unsigned int firstGap = aln.qStartPos;
        while (firstGap) // need to make as much 127-long deletion blocks as needed
        {
            unsigned int gapToWrite = -std::min((unsigned int) (127), firstGap);
            firstBlock.write(reinterpret_cast<const char *>(&zero),
                             sizeof(char)); // no match in this block,
            firstBlock.write(reinterpret_cast<const char *>(&gapToWrite), sizeof(char)); // only gaps
            firstGap -= -gapToWrite;
            nbOfBlocks++;
        }


        unsigned short int startPos = aln.dbStartPos + 1; // Starts at 1 in ca3m format
        unsigned int targetKey = i ? referenceDBr->dbBKeyMap(aln.dbKey) : referenceDBr->dbAKeyMap(aln.dbKey);
        unsigned int targetIdentifier = referenceDBr->getId(targetKey);
        out.write(reinterpret_cast<const char *>(&targetIdentifier), sizeof(unsigned int));
        out.write(reinterpret_cast<const char *>(&startPos), sizeof(unsigned short int));
        out.write(reinterpret_cast<const char *>(&nbOfBlocks), sizeof(unsigned short int));
        out << firstBlock.str() << blocksDescription.str();

        totalNumberOfBlocks += nbOfBlocks;
    }

    return out.str();
}
