/*
 * extractA3M adapted from HH-suite 3.0: a3m_compress.h
 * Original Author: meiermark
 * Licensed under GPLv3
 */

#include "CompressedA3M.h"
#include "DBConcat.h"
#include <sstream>

void readU16(const char **ptr, uint16_t &result) {
    unsigned char array[2];

    array[0] = (unsigned char) (**ptr);
    (*ptr)++;
    array[1] = (unsigned char) (**ptr);
    (*ptr)++;

    result = array[0] | (array[1] << 8);
}


void readU32(const char **ptr, uint32_t &result) {
    unsigned char array[4];

    array[0] = (unsigned char) (**ptr);
    (*ptr)++;
    array[1] = (unsigned char) (**ptr);
    (*ptr)++;
    array[2] = (unsigned char) (**ptr);
    (*ptr)++;
    array[3] = (unsigned char) (**ptr);
    (*ptr)++;

    result = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);
}

std::string CompressedA3M::extractA3M(const char *data, size_t data_size,
                                      DBReader<unsigned int>& sequenceReader,
                                      DBReader<unsigned int>& headerReader) {
    std::ostringstream output;

    //read stuff till compressed part
    char last_char = '\0';
    size_t index = 0;
    size_t consensus_length = 0;
    char inConsensus = 0;

    //handle commentary line
    if ((*data) == '#') {
        while ((*data) != '\n') {
            output.put((*data));
            last_char = (*data);
            data++;
            index++;
        }

        output.put('\n');
        last_char = '\n';
        data++;
        index++;
    }

    while (!(last_char == '\n' && (*data) == ';') && index < data_size) {
        if ((*data) == '\n') {
            inConsensus++;
        }
        else if (inConsensus == 1) {
            consensus_length++;
        }

        output.put((*data));
        last_char = (*data);
        data++;
        index++;
    }

    //get past ';'
    data++;
    index++;

    while (index < data_size - 1) {
        unsigned int entry_index;
        unsigned short int nr_blocks;
        unsigned short int start_pos;

        readU32(&data, entry_index);
        index += 4;

        std::string sequence = sequenceReader.getData(entry_index);
        std::string header = headerReader.getData(entry_index);

        // make sure we always have a valid fasta prefix
        if (header[0] != '>') {
            output.put('>');
        }

        output.write(header.c_str(), header.length() - 1);
        output.put('\n');

        readU16(&data, start_pos);
        index += 2;

        readU16(&data, nr_blocks);
        index += 2;

        size_t actual_pos = start_pos;
        size_t alignment_length = 0;
        for (unsigned short int block_index = 0; block_index < nr_blocks; block_index++) {
            unsigned char nr_matches = (unsigned char) (*data);
            data++;
            index++;

            for (int i = 0; i < nr_matches; i++) {
                output.put(sequence[actual_pos - 1]);
                actual_pos++;
                alignment_length++;
            }

            char nr_insertions_deletions = (*data);
            data++;
            index++;

            if (nr_insertions_deletions > 0) {
                for (int i = 0; i < nr_insertions_deletions; i++) {
                    output.put(tolower(sequence[actual_pos - 1]));
                    actual_pos++;
                }
            }
            else {
                for (int i = 0; i < -nr_insertions_deletions; i++) {
                    output.put('-');
                    alignment_length++;
                }
            }
        }

        while (alignment_length < consensus_length) {
            output.put('-');
            alignment_length++;
        }

        output.put('\n');
    }

    return output.str();
}


std::string CompressedA3M::fromAlignmentResult(const std::vector<Matcher::result_t>& alignment, DBConcat& referenceDBr) {
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
        unsigned int targetKey = i ? referenceDBr.dbBKeyMap(aln.dbKey) : referenceDBr.dbAKeyMap(aln.dbKey);
        unsigned int targetIdentifier = referenceDBr.getId(targetKey);
        out.write(reinterpret_cast<const char *>(&targetIdentifier), sizeof(unsigned int));
        out.write(reinterpret_cast<const char *>(&startPos), sizeof(unsigned short int));
        out.write(reinterpret_cast<const char *>(&nbOfBlocks), sizeof(unsigned short int));
        out << firstBlock.str() << blocksDescription.str();

        totalNumberOfBlocks += nbOfBlocks;
    }

    return out.str();
}
