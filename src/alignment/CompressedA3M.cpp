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
                                      DBReader<unsigned int>& headerReader, int thread_idx) {
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
        uint32_t entry_index;
        uint16_t nr_blocks;
        uint16_t start_pos;

        readU32(&data, entry_index);
        index += 4;

        std::string sequence = sequenceReader.getData(entry_index, thread_idx);
        std::string header = headerReader.getData(entry_index, thread_idx);

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

void CompressedA3M::extractMatcherResults(unsigned int &key, std::vector<Matcher::result_t> &results,
        const char *data, size_t dataSize, DBReader<unsigned int> &sequenceReader, bool skipFirst) {
    //read stuff till compressed part
    char lastChar = '\0';
    size_t index = 0;

    //handle commentary line
    if ((*data) == '#') {
        while ((*data) != '\n') {
            data++;
            index++;
        }

        lastChar = '\n';
        data++;
        index++;
    }

    char inConsensus = 0;
    size_t consensusLength = 0;
    while (!(lastChar == '\n' && (*data) == ';') && index < dataSize) {
        if ((*data) == '\n') {
            inConsensus++;
        }
        else if (inConsensus == 1) {
            consensusLength++;
        }

        lastChar = (*data);
        data++;
        index++;
    }

    std::string backtrace;
    backtrace.reserve(consensusLength);

    //get past ';'
    data++;
    index++;

    size_t qLen = 0;
    bool isFirst = true;
    while (index < dataSize - 1) {
        Matcher::result_t match;
        match.seqId      = 0;
        match.score      = 0;
        match.eval       = 0;

        uint32_t entryIndex;
        readU32(&data, entryIndex);
        index += 4;

        match.dbKey = sequenceReader.getDbKey(entryIndex);
        if (isFirst) {
            key = match.dbKey;
            qLen = sequenceReader.getSeqLen(entryIndex);
            match.qLen = match.dbLen = qLen;
        } else {
            match.qLen = qLen;
            match.dbLen = sequenceReader.getSeqLen(entryIndex);
        }

        unsigned short int startPos;
        readU16(&data, startPos);
        index += 2;

        match.dbStartPos = startPos - 1;

        unsigned short int nrBlocks;
        readU16(&data, nrBlocks);
        index += 2;

        if (skipFirst && isFirst) {
            data += 2 * nrBlocks;
            index += 2 * nrBlocks;
            isFirst = false;
            continue;
        }

        match.qStartPos = 0;
        unsigned int qAlnLength = 0;
        unsigned int dbAlnLength = 0;
        bool firstBlockM = false;
        for (unsigned short int blockIdx = 0; blockIdx < nrBlocks; blockIdx++) {
            unsigned char matchCount = (unsigned char) (*data);
            data++;
            index++;

            qAlnLength += matchCount;
            dbAlnLength += matchCount;
            backtrace.append(matchCount, 'M');

            if (matchCount != 0) {
                firstBlockM = true;
            }

            signed char inDelCount = (*data);
            data++;
            index++;

            if (firstBlockM == false) {
                match.qStartPos -= inDelCount;
            } else {
                if (inDelCount > 0) {
                    backtrace.append(inDelCount, 'D');
                    qAlnLength += inDelCount;
                } else if (inDelCount < 0) {
                    backtrace.append(-inDelCount, 'I');
                    dbAlnLength -= inDelCount;
                }
            }
        }

        match.backtrace = backtrace;
        match.qEndPos = match.qStartPos + dbAlnLength - 1;
        match.dbEndPos = match.dbStartPos + qAlnLength - 1;
        results.emplace_back(match);
        backtrace.clear();
    }
}


//#define CA3M_DEBUG
void CompressedA3M::hitToBuffer(unsigned int targetId, const Matcher::result_t& hit, std::string& buffer) {
#ifndef CA3M_DEBUG
    buffer.append(reinterpret_cast<const char *>(&targetId), sizeof(unsigned int));
#else
    buffer.append(SSTR((int)targetId));
    buffer.append(1, '\t');
#endif
    // Starts at 1 in ca3m format
    unsigned short int startPos = hit.dbStartPos + 1;
#ifndef CA3M_DEBUG
    buffer.append(reinterpret_cast<const char *>(&startPos), sizeof(unsigned short));
#else
    buffer.append(SSTR((int)startPos));
    buffer.append(1, '\t');
#endif
    unsigned int nbOfBlocks = 0;
    size_t countPos = buffer.size();
#ifndef CA3M_DEBUG
    buffer.append(sizeof(unsigned short), '\0');
#endif
    // add the deletion at the begining of the sequence (local alignment)
    unsigned int firstGap = hit.qStartPos;
    // need to make as many 127-long deletion blocks as needed
    while (firstGap) {
        unsigned int gapToWrite = -std::min((unsigned int) (127), firstGap);
#ifndef CA3M_DEBUG
        // no match in this block
        buffer.append(sizeof(char), '\0');
        // only gaps
        buffer.append(reinterpret_cast<const char *>(&gapToWrite), sizeof(char));
#else
        buffer.append(SSTR((int)0));
        buffer.append(1, '\t');
        buffer.append(SSTR((int)gapToWrite));
        buffer.append(1, '\t');
#endif

        firstGap -= -gapToWrite;
        nbOfBlocks++;
    }

    // detect the blocks
    for (size_t btIndex = 0; btIndex < hit.backtrace.size();) {
        // seek the next insertion or deletion
        unsigned char matchLen = 0;
        while (btIndex < hit.backtrace.size() && hit.backtrace[btIndex] == 'M' && matchLen < 255) {
            btIndex++;
            matchLen++;
        }

#ifndef CA3M_DEBUG
        buffer.append(reinterpret_cast<const char *>(&matchLen), sizeof(unsigned char));
#else
        buffer.append(SSTR((int)matchLen));
        buffer.append(1, '\t');
#endif
        // end of block because an I or D was found
        char inOrDel = 0;
        if (btIndex < hit.backtrace.size() && hit.backtrace[btIndex] != 'M') {
            // store whether it is I or D
            inOrDel = hit.backtrace[btIndex];
        }

        // seek the next match
        unsigned char indelLen = 0;
        while (btIndex < hit.backtrace.size() && hit.backtrace[btIndex] == inOrDel && indelLen < 127) {
            btIndex++;
            indelLen++;
        }

        // deletions must be counted backwards
        // deletion states are insertion states from the perspective of the aligned sequence
        if (indelLen && inOrDel == 'I') {
            indelLen *= -1;
        }
#ifndef CA3M_DEBUG
        buffer.append(sizeof(char), indelLen);
#else
        buffer.append(SSTR((int)indelLen));
        buffer.append(1, '\t');
#endif
        nbOfBlocks++;
    }
#ifndef CA3M_DEBUG
    buffer.replace(countPos, sizeof(unsigned short), reinterpret_cast<const char *>(&nbOfBlocks), sizeof(unsigned short));
#else
    buffer.append(SSTR((int)nbOfBlocks));
    buffer.append(1, '\n');
#endif
}
