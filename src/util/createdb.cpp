/*
 * createdb
 * written by Martin Steinegger <martin.steinegger@snu.ac.kr>.
 * modified by Maria Hauser <mhauser@genzentrum.lmu.de> (splitting into sequences/headers databases)
 * modified by Milot Mirdita <milot@mirdita.de>
 */
#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "KSeqWrapper.h"
#include "itoa.h"
#include "FastSort.h"
#include "Masker.h"
#include "MMseqsNode.h"

#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

// Sort the data file in-place using your index array
int sortWithIndex(const char *dataFileSeq,
                  const char *indexFileSeq,
                  const char *dataFileHeader,
                  const char *indexFileHeader)
{
    DBReader<DBKeyType> reader(dataFileSeq, indexFileSeq, 1, DBReader<DBKeyType>::USE_INDEX);
    reader.open(DBReader<DBKeyType>::HARDNOSORT);
    DBReader<DBKeyType>::Index *index = reader.getIndex();
    struct stat st;
    if (stat(dataFileSeq, &st) != 0) {
        Debug(Debug::ERROR) << "stat failed: " << dataFileSeq << "\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t seqfileSize = st.st_size;

    FILE *fin = FileUtil::openFileOrDie(dataFileSeq, "rb", true);
    char * buf = new char [seqfileSize];
    size_t got = fread(buf, 1, seqfileSize, fin);
    if (got != seqfileSize) {
        Debug(Debug::ERROR) << "short read: " << dataFileSeq << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(fin);
    // needed to keep the information in what line the id was originally
    for (size_t i = 0; i < reader.getSize(); i++) {
        index[i].id = i;
    }

    SORT_PARALLEL(index, index + reader.getSize(), DBReader<DBKeyType>::Index::compareByLength);

    FILE *seqOut = FileUtil::openFileOrDie(dataFileSeq, "wb", true);
    setvbuf(seqOut, NULL, _IOFBF, 1024*1024*50);
    size_t offset = 0;
    for (size_t i = 0; i < reader.getSize(); i++) {
        size_t written = fwrite(buf + index[i].offset, 1, index[i].length, seqOut);
        if (written != static_cast<size_t>(index[i].length)) {
            Debug(Debug::ERROR) << "Can not write to data file " << dataFileSeq << "\n";
            EXIT(EXIT_FAILURE);
        }
        index[i].offset = offset;
        offset += written;
    }
    fclose(seqOut);

    if (stat(dataFileHeader, &st) != 0) {
        Debug(Debug::ERROR) << "stat failed: " << dataFileHeader << "\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t headFileSize = st.st_size;
    if(headFileSize > seqfileSize){
        delete [] buf;
        buf = new char [headFileSize];
    }
    fin = FileUtil::openFileOrDie(dataFileHeader, "rb", true);
    got = fread(buf, 1, headFileSize, fin);
    if (got != headFileSize) {
        Debug(Debug::ERROR) << "short read: " << dataFileHeader << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(fin);

    DBReader<DBKeyType> header(dataFileHeader, indexFileHeader, 1, DBReader<DBKeyType>::USE_INDEX);
    header.open(DBReader<DBKeyType>::HARDNOSORT);
    DBReader<DBKeyType>::Index *headerIndex = header.getIndex();
    FILE *headerout = FileUtil::openFileOrDie(dataFileHeader, "wb", true);
    setvbuf(headerout, NULL, _IOFBF, 1024*1024*50);
    offset = 0;
    for (size_t i = 0; i < header.getSize(); i++) {
        size_t sortedId = static_cast<size_t>(index[i].id);
        size_t written = fwrite(buf + headerIndex[sortedId].offset, 1, headerIndex[sortedId].length, headerout);
        // reconstruct old id
        index[i].id = headerIndex[sortedId].id;
        if (written != static_cast<size_t>(headerIndex[sortedId].length)) {
            Debug(Debug::ERROR) << "Can not write to data file " << dataFileHeader << "\n";
            EXIT(EXIT_FAILURE);
        }
        headerIndex[sortedId].offset = offset;
        offset += written;
    }
    fclose(headerout);
    delete [] buf;

    SORT_PARALLEL(index, index + reader.getSize(), DBReader<DBKeyType>::Index::compareByOffset);
    {
        std::string tmpIndex = std::string(indexFileSeq) + ".tmp";
        FILE *indexout = FileUtil::openFileOrDie(tmpIndex.c_str(), "wb", false);
        DBWriter::writeIndex(indexout, reader.getSize(), index);
        fclose(indexout);
        FileUtil::move(tmpIndex.c_str(), indexFileSeq);

        std::string tmpHeaderIndex = std::string(indexFileHeader) + ".tmp";
        FILE *headerIndexOut = FileUtil::openFileOrDie(tmpHeaderIndex.c_str(), "wb", false);
        DBWriter::writeIndex(headerIndexOut, header.getSize(), headerIndex);
        fclose(headerIndexOut);
        FileUtil::move(tmpHeaderIndex.c_str(), indexFileHeader);
    }

    reader.close();
    header.close();
    return 0;
}

int mergeSequentialByJointIndex(
        char ** dataFiles,
        char ** indexFiles,
        char ** dataFilesHeader,
        char ** indexFilesHeader,
        char * outDataFile,
        char * outIndexFile,
        char * outHeaderDataFile,
        char * outHeaderIndexFile,
        const char * outLookupFile,
        std::vector<unsigned int>* sourceLookup,
        size_t totalEntries,
        size_t shuffleSplits
) {
    struct JointEntry {
        unsigned int fileIdx;
        DBKeyType id;
        unsigned int length;
        JointEntry(unsigned int fileIdx, DBKeyType id, unsigned int length)
            : fileIdx(fileIdx), id(id), length(length) {}

        bool operator<(JointEntry const &o) const {
            if (length != o.length){
                return length < o.length;
            }
            return id < o.id;
        }
    };

    std::vector<JointEntry> joint;
    joint.reserve(totalEntries);
    size_t maxLen = 0;
    for (size_t i = 0; i < shuffleSplits; i++) {
        DBReader<DBKeyType> reader(
                dataFiles[i],
                indexFiles[i],
                1,
                DBReader<DBKeyType>::USE_INDEX
        );
        reader.open(DBReader<DBKeyType>::HARDNOSORT);
        DBReader<DBKeyType>::Index* index = reader.getIndex();
        for(size_t j = 0; j < reader.getSize(); j++){
            joint.emplace_back((unsigned int)i, index[j].id, index[j].length);
            maxLen = std::max(maxLen, static_cast<size_t>(index[j].length));
        }
        reader.close();
    }

    SORT_PARALLEL(joint.begin(), joint.end());

    // 4) Open each data file once (no fseek later)
    std::vector<FILE*> inFileSeq(shuffleSplits);
    std::vector<FILE*> inFileHeader(shuffleSplits);

    for (size_t i = 0; i < shuffleSplits; i++) {
        inFileSeq[i] = FileUtil::openFileOrDie(
            dataFiles[i], "rb", true
        );
        inFileHeader[i] = FileUtil::openFileOrDie(
            dataFilesHeader[i], "rb", true
        );
    }

    FILE *fout = FileUtil::openAndDelete(outDataFile, "wb");
    setvbuf(fout, NULL, _IOFBF, 1024*1024*50);
    FILE *idxOut = FileUtil::openAndDelete(outIndexFile, "wb");
    setvbuf(idxOut, NULL, _IOFBF, 1024*1024*50);
    FILE *foutHeader  = FileUtil::openAndDelete(outHeaderDataFile, "wb");
    setvbuf(foutHeader, NULL, _IOFBF, 1024*1024*50);
    FILE *foutLookup  = FileUtil::openAndDelete(outLookupFile, "wb");
    setvbuf(foutLookup, NULL, _IOFBF, 1024*1024*50);
    FILE *idxOutHeader = FileUtil::openAndDelete(outHeaderIndexFile, "wb");
    setvbuf(idxOutHeader, NULL, _IOFBF, 1024*1024*50);

    struct FileBuffer {
        std::vector<char> buffer;
        size_t pos, end;
        FILE *file;
        FileBuffer(FILE* f) : buffer(1024 * 1024), pos(0), end(0), file(f) {}

        bool refill() {
            pos = 0;
            // read 1 MB
            end = fread(buffer.data(), 1, 1024 * 1024, file);
            return end != 0;
        }

        bool getChar(char &ch) {
            if (pos >= end && !refill()) { return false; }
            ch = buffer[pos++];
            return true;
        }
    };

    // Create buffers for each header file
    std::vector<FileBuffer> headerBuffers;
    std::vector<char> writeHeaderBuf;
    headerBuffers.reserve(shuffleSplits);
    for (size_t i = 0; i < shuffleSplits; i++) {
        headerBuffers.emplace_back(inFileHeader[i]);
    }
    size_t mergedOffset = 0;
    size_t mergedOffsetHeader = 0;
    std::vector<char> scratch(maxLen);
    DBReader<DBKeyType>::Index entry;
    DBReader<DBKeyType>::LookupEntry lookupEntry;

    char indexBuffer[1024];
    std::string lookupBuffer;
    const int ALIGN = 4;
    const char pad_buffer[4] = {20, 20, 20, 20}; // pre-filled buffer
    for(size_t i = 0; i < joint.size(); i++){
        JointEntry qe = joint[i];
        size_t read = fread(scratch.data(), 1, qe.length, inFileSeq[qe.fileIdx]);
        if (UNLIKELY(read != static_cast<size_t>(qe.length))) {
            Debug(Debug::ERROR) << "Can not read to data file " << dataFiles[qe.fileIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }
        size_t written = fwrite(scratch.data(), 1, qe.length, fout);
        const size_t sequencepadding = (qe.length % ALIGN == 0) ? 0 : ALIGN - qe.length % ALIGN;
        written +=  fwrite(pad_buffer, 1, sequencepadding, fout);
        if (UNLIKELY(written != qe.length + sequencepadding)) {
            Debug(Debug::ERROR) << "Can not write to data file " << outDataFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        writeHeaderBuf.clear();
        char ch;
        do {
            if (!headerBuffers[qe.fileIdx].getChar(ch)) {
                Debug(Debug::ERROR) << "Unexpected EOF in header file " <<  qe.fileIdx << "\n";
                EXIT(EXIT_FAILURE);
            }
            writeHeaderBuf.push_back(ch);
        } while (ch != '\0');

        lookupEntry.id = i;
        lookupEntry.entryName = Util::parseFastaHeader(writeHeaderBuf.data());
        if (lookupEntry.entryName.empty()) {
            Debug(Debug::WARNING) << "Cannot extract identifier from entry " << lookupEntry.id  << "\n";
        }
        lookupEntry.fileNumber = sourceLookup[qe.fileIdx][(qe.id - qe.fileIdx) / 32];
        lookupBuffer.clear();
        DBReader<DBKeyType>::lookupEntryToBuffer(lookupBuffer, lookupEntry);
        written = fwrite(lookupBuffer.data(), 1, lookupBuffer.size(), foutLookup);
        if (UNLIKELY(written != lookupBuffer.size())) {
            Debug(Debug::ERROR) << "Can not write to lookup file " << outLookupFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        written = fwrite(writeHeaderBuf.data(), 1, writeHeaderBuf.size(), foutHeader);
        if (UNLIKELY(written != writeHeaderBuf.size())) {
            Debug(Debug::ERROR) << "Can not write to header file " << outHeaderDataFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        entry.offset = mergedOffset;
        // + 2 is needed for newline and null character
        entry.length = qe.length + 2;
        entry.id = i;
        DBWriter::writeIndexEntryToFile(idxOut, indexBuffer, entry);
        entry.length = writeHeaderBuf.size();
        entry.offset = mergedOffsetHeader;
        DBWriter::writeIndexEntryToFile(idxOutHeader, indexBuffer, entry);
        mergedOffset += qe.length + sequencepadding;
        mergedOffsetHeader += writeHeaderBuf.size();
    }

    fclose(fout);
    fclose(idxOut);
    fclose(foutHeader);
    fclose(foutLookup);
    fclose(idxOutHeader);
    for (size_t i = 0; i < inFileSeq.size(); i++) {
        fclose(inFileSeq[i]);
        fclose(inFileHeader[i]);
        FileUtil::remove(dataFiles[i]);
        FileUtil::remove(indexFiles[i]);
        FileUtil::remove(dataFilesHeader[i]);
        FileUtil::remove(indexFilesHeader[i]);
    }
    return 0;
}

void processSeqBatch(Parameters & par, DBWriter &seqWriter, DBWriter &hdrWriter, BaseMatrix *subMat, int querySeqType,
                     Masker ** masker, Sequence ** seqs, size_t currId,
                     std::vector<std::pair<std::vector<char>, std::string>> &entries, const size_t entriesSize,
                     unsigned int shuffleSplits){
    if(masker[0] == NULL){
        for(int i = 0; i < par.threads; i++){
            masker[i] = new Masker(*subMat);
            seqs[i] = new Sequence(par.maxSeqLen, querySeqType, subMat, 0, false, false);
        }
    }

#pragma omp parallel num_threads(par.threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < entriesSize; i++) {
            seqs[thread_idx]->mapSequence(currId + i, currId + i, entries[i].first.data(), entries[i].first.size());
            const unsigned char *numSequence = seqs[thread_idx]->numSequence;
            std::copy_n(numSequence, seqs[thread_idx]->L, entries[i].first.begin());
            masker[thread_idx]->maskSequence(*seqs[thread_idx], par.maskMode, par.maskProb, par.maskLowerCaseMode,
                                             par.maskNrepeats);

            for (int j = 0; j < seqs[thread_idx]->L; j++) {
                entries[i].first[j] = (numSequence[j] == masker[thread_idx]->maskLetterNum) ? entries[i].first[j] + 32
                                                                                            : entries[i].first[j];
            }
        }
    }

    for (size_t i = 0; i < entriesSize; i++) {
        size_t id = currId + i;
        size_t splitIdx = id % shuffleSplits;
        seqWriter.writeData(entries[i].first.data(), entries[i].first.size(), currId + i, splitIdx, false);
        hdrWriter.writeData(entries[i].second.c_str(), entries[i].second.length(), currId + i, splitIdx);
    }
}

struct LengthBin {
    size_t length;
    size_t sequenceCount;
    size_t headerSize;
    size_t startKey;
    size_t offset;
    size_t headerOffset;
    size_t file;
};

struct InputPiece {
    size_t file;
    size_t start;
    size_t end;
};

std::vector<LengthBin> makeLengthBins(const std::vector<size_t> &sequenceCount, const std::vector<size_t> &headerSize) {
    std::vector<LengthBin> histogram;
    for (size_t length = sequenceCount.size(); length-- > 0; ) {
        if (sequenceCount[length] > 0) {
            LengthBin bin = {length, sequenceCount[length], headerSize[length], 0, 0, 0, 0};
            histogram.push_back(bin);
        }
    }
    return histogram;
}

std::vector<LengthBin> mergeLengthHistograms(const std::vector<std::vector<LengthBin>> &histograms) {
    size_t maxLength = 0;
    for (size_t i = 0; i < histograms.size(); i++) {
        if (histograms[i].empty() == false) {
            maxLength = std::max(maxLength, histograms[i][0].length);
        }
    }
    std::vector<size_t> sequenceCount(maxLength + 1, 0);
    std::vector<size_t> headerSize(maxLength + 1, 0);
    for (size_t i = 0; i < histograms.size(); i++) {
        for (size_t j = 0; j < histograms[i].size(); j++) {
            sequenceCount[histograms[i][j].length] += histograms[i][j].sequenceCount;
            headerSize[histograms[i][j].length] += histograms[i][j].headerSize;
        }
    }
    return makeLengthBins(sequenceCount, headerSize);
}

std::vector<LengthBin> getLengthToBin(const std::vector<LengthBin> &histogram) {
    std::vector<LengthBin> lengthToBin(histogram.empty() ? 0 : histogram[0].length + 1);
    for (size_t i = 0; i < histogram.size(); i++) {
        lengthToBin[histogram[i].length] = histogram[i];
    }
    return lengthToBin;
}

void splitStarts(const std::vector<LengthBin> &starts, std::vector<std::vector<LengthBin>> &splits) {
    std::vector<LengthBin> lengthToStart = getLengthToBin(starts);
    for (size_t i = 0; i < splits.size(); i++) {
        for (size_t j = 0; j < splits[i].size(); j++) {
            LengthBin &bin = splits[i][j];
            LengthBin &start = lengthToStart[bin.length];
            bin.startKey = start.startKey;
            bin.offset = start.offset;
            bin.headerOffset = start.headerOffset;
            bin.file = start.file;
            start.startKey += bin.sequenceCount;
            start.offset += bin.sequenceCount * (bin.length + 2);
            start.headerOffset += bin.headerSize;
        }
    }
}

void flushEntryBuffer(const std::vector<int> &dataFileFds, const std::string &dataFile, const char *buffer, size_t bufferPos, char *sortedBuffer, std::vector<LengthBin> &lengthToBin) {
    const size_t lengthCount = lengthToBin.size();
    std::vector<size_t> runStart(2 * lengthCount + 1, 0);
    for (size_t pos = 0; pos < bufferPos; ) {
        const size_t seqSize = strlen(buffer + pos) + 1;
        const size_t hdrSize = strlen(buffer + pos + seqSize) + 1;
        runStart[seqSize - 1] += seqSize;
        runStart[lengthCount + seqSize - 1] += hdrSize;
        pos += seqSize + hdrSize;
    }
    for (size_t i = 1; i < runStart.size(); i++) {
        runStart[i] += runStart[i - 1];
    }
    std::vector<size_t> runPos(runStart.begin(), runStart.end() - 1);
    for (size_t pos = 0; pos < bufferPos; ) {
        const size_t seqSize = strlen(buffer + pos) + 1;
        const size_t hdrSize = strlen(buffer + pos + seqSize) + 1;
        const size_t length = seqSize - 2;
        memcpy(sortedBuffer + runPos[length], buffer + pos, seqSize);
        runPos[length] += seqSize;
        memcpy(sortedBuffer + runPos[lengthCount + length], buffer + pos + seqSize, hdrSize);
        runPos[lengthCount + length] += hdrSize;
        pos += seqSize + hdrSize;
    }
    for (size_t i = 0; i < 2 * lengthCount; i++) {
        const bool isHeader = i >= lengthCount;
        LengthBin &bin = lengthToBin[isHeader ? i - lengthCount : i];
        size_t &offset = isHeader ? bin.headerOffset : bin.offset;
        const size_t size = runStart[i + 1] - runStart[i];
        for (size_t written = 0; written < size; ) {
            const ssize_t count = pwrite(dataFileFds[isHeader ? MMseqsNode::count + bin.file : bin.file], sortedBuffer + runStart[i] + written, size - written, offset + written);
            if (count <= 0) {
                Debug(Debug::ERROR) << "Can not write to data file " << dataFile << (isHeader ? "_h." : ".") << bin.file << "\n";
                EXIT(EXIT_FAILURE);
            }
            written += count;
        }
        offset += size;
    }
}

template <typename T>
std::vector<T> readLinfinityFile(const std::string &file) {
    const std::string linfinityData = MMseqsNode::waitForFile(file);
    char *linfinityDataChar = (char *) linfinityData.c_str();
    const size_t linfinityDataSize = linfinityData.size();
    const size_t columns = sizeof(T) / sizeof(size_t);
    std::vector<T> entries;
    size_t currPos = 0;
    const char *cols[columns];
    while (currPos < linfinityDataSize) {
        const size_t colCount = Util::getWordsOfLine(linfinityDataChar, cols, columns);
        size_t values[columns] = {};
        for (size_t i = 0; i < colCount; i++) {
            values[i] = Util::fast_atoi<size_t>(cols[i]);
        }
        T entry;
        memcpy(&entry, values, sizeof(T));
        entries.push_back(entry);
        linfinityDataChar = Util::skipLine(linfinityDataChar);
        currPos = linfinityDataChar - linfinityData.c_str();
    }
    return entries;
}

int createdbLinfinity(Parameters &par, const std::vector<std::string> &filenames, const std::string &dataFile) {
    const std::string planPrefix = dataFile + ".splitplan.";
    if (MMseqsNode::isMaster()) {
        std::vector<size_t> fileSizes(filenames.size());
        size_t totalSize = 0;
        for (size_t i = 0; i < filenames.size(); i++) {
            fileSizes[i] = FileUtil::getFileSize(filenames[i]);
            if (fileSizes[i] == SIZE_MAX) {
                Debug(Debug::ERROR) << "File " << filenames[i] << " does not exist\n";
                EXIT(EXIT_FAILURE);
            }
            totalSize += fileSizes[i];
        }
        const size_t pieceSize = std::max<size_t>(1, totalSize / (MMseqsNode::count * par.threads));
        std::vector<std::vector<size_t>> fileCuts(filenames.size(), std::vector<size_t>(1, 0));
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (size_t i = 0; i < filenames.size(); i++) {
            std::vector<size_t> &cuts = fileCuts[i];
            if (fileSizes[i] > pieceSize && Util::endsWith(".gz", filenames[i]) == false && Util::endsWith(".bz2", filenames[i]) == false) {
                KSeqWrapper *kseq = KSeqFactory(filenames[i].c_str());
                const bool isFasta = kseq->ReadEntry() && kseq->entry.qual.l == 0;
                delete kseq;
                if (isFasta && Util::endsWith(".zst", filenames[i]) == true) {
                    FILE *handle = FileUtil::openFileOrDie(filenames[i].c_str(), "r", true);
                    setvbuf(handle, NULL, _IONBF, 0);
                    uint32_t header[3];
                    size_t frame = 0;
                    while (frame < fileSizes[i] && fseek(handle, frame, SEEK_SET) == 0 && fread(header, sizeof(uint32_t), 3, handle) == 3
                           && header[0] == ZSTD_MAGIC_SKIPPABLE_START && header[1] == 4) {
                        if (frame - cuts.back() >= pieceSize) {
                            cuts.push_back(frame);
                        }
                        frame += sizeof(header) + header[2];
                    }
                    if (fclose(handle) != 0) {
                        Debug(Debug::ERROR) << "Cannot close file " << filenames[i] << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    if (frame != fileSizes[i]) {
                        cuts.resize(1);
                    }
                } else if (isFasta) {
                    for (size_t cut = pieceSize; cut < fileSizes[i]; cut += pieceSize) {
                        cuts.push_back(cut);
                    }
                }
            }
            cuts.push_back(fileSizes[i]);
        }
        std::vector<InputPiece> pieces;
        std::vector<size_t> pieceOffset(1, 0);
        for (size_t i = 0; i < filenames.size(); i++) {
            for (size_t j = 0; j + 1 < fileCuts[i].size(); j++) {
                InputPiece piece = {i, fileCuts[i][j], fileCuts[i][j + 1]};
                pieces.push_back(piece);
                pieceOffset.push_back(pieceOffset.back() + (piece.end - piece.start));
            }
        }
        const size_t pieceCount = pieces.size();
        std::vector<size_t> pieceBoundaries(MMseqsNode::count + 1, pieceCount);
        pieceBoundaries[0] = 0;
        for (unsigned int node = 0; node + 1 < MMseqsNode::count; node++) {
            const size_t targetSize = totalSize / MMseqsNode::count * (node + 1);
            const size_t maxPieceEnd = (pieceCount >= MMseqsNode::count) ? pieceCount - (MMseqsNode::count - node - 1) : pieceCount;
            size_t nearestEnd = std::lower_bound(pieceOffset.begin(), pieceOffset.end(), targetSize) - pieceOffset.begin();
            if (nearestEnd > 0 && targetSize - pieceOffset[nearestEnd - 1] <= pieceOffset[nearestEnd] - targetSize) {
                nearestEnd--;
            }
            pieceBoundaries[node + 1] = std::min(std::max(nearestEnd, pieceBoundaries[node] + 1), maxPieceEnd);
        }
        for (unsigned int node = 0; node < MMseqsNode::count; node++) {
            std::string planData;
            std::string nodeFileNames;
            for (size_t i = pieceBoundaries[node]; i < pieceBoundaries[node + 1]; i++) {
                planData += SSTR(pieces[i].file) + "\t" + SSTR(pieces[i].start) + "\t" + SSTR(pieces[i].end) + "\n";
                if (i == pieceBoundaries[node] || pieces[i].file != pieces[i - 1].file) {
                    nodeFileNames.append(1, ' ');
                    nodeFileNames.append(filenames[pieces[i].file]);
                }
            }
            MMseqsNode::publishFile(planPrefix + SSTR(node), planData);
            Debug(Debug::INFO) << "Node " << node << " reads" << nodeFileNames << "\n";
        }
    }
    const std::string planFile = planPrefix + SSTR(MMseqsNode::rank);
    const std::vector<InputPiece> nodePieces = readLinfinityFile<InputPiece>(planFile);
    FileUtil::remove(planFile.c_str());

    Debug(Debug::INFO) << "Converting sequences\n";
    std::vector<std::vector<LengthBin>> pieceHistograms(nodePieces.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
    for (size_t i = 0; i < nodePieces.size(); i++) {
        std::vector<size_t> sequenceCount;
        std::vector<size_t> headerSize;
        KSeqWrapper *kseq = KSeqFactory(filenames[nodePieces[i].file].c_str(), nodePieces[i].start, nodePieces[i].end);
        while (kseq->ReadEntry()) {
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            const size_t nameLength = strlen(e.name.s);
            if (nameLength == 0) {
                Debug(Debug::ERROR) << "Fasta entry in " << filenames[nodePieces[i].file] << " is invalid\n";
                EXIT(EXIT_FAILURE);
            }
            const size_t length = e.sequence.l;
            if (length >= sequenceCount.size()) {
                sequenceCount.resize(length + 1, 0);
                headerSize.resize(length + 1, 0);
            }
            sequenceCount[length]++;
            headerSize[length] += nameLength + 2;
        }
        delete kseq;
        pieceHistograms[i] = makeLengthBins(sequenceCount, headerSize);
    }

    const std::vector<LengthBin> localNodeHistogram = mergeLengthHistograms(pieceHistograms);

    const std::string histogramPrefix = dataFile + ".histogram.";
    const std::string startsPrefix = dataFile + ".starts.";
    std::string histogramData;
    for (size_t i = 0; i < localNodeHistogram.size(); i++) {
        histogramData += SSTR(localNodeHistogram[i].length) + "\t" + SSTR(localNodeHistogram[i].sequenceCount) + "\t" + SSTR(localNodeHistogram[i].headerSize) + "\n";
    }
    MMseqsNode::publishFile(histogramPrefix + SSTR(MMseqsNode::rank), histogramData);
    if (MMseqsNode::isMaster()) {
        std::vector<std::vector<LengthBin>> nodeHistograms(MMseqsNode::count);
        for (unsigned int node = 0; node < MMseqsNode::count; node++) {
            const std::string histogramFile = histogramPrefix + SSTR(node);
            nodeHistograms[node] = readLinfinityFile<LengthBin>(histogramFile);
            FileUtil::remove(histogramFile.c_str());
        }

        std::vector<LengthBin> globalHistogram = mergeLengthHistograms(nodeHistograms);
        if (globalHistogram.empty() == false && globalHistogram.back().length == 0) {
            Debug(Debug::WARNING) << globalHistogram.back().sequenceCount << " sequences without residues are skipped\n";
            globalHistogram.pop_back();
            for (unsigned int node = 0; node < MMseqsNode::count; node++) {
                if (nodeHistograms[node].empty() == false && nodeHistograms[node].back().length == 0) {
                    nodeHistograms[node].pop_back();
                }
            }
        }
        if (globalHistogram.empty()) {
            Debug(Debug::ERROR) << "The input files have no entry: ";
            for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
                Debug(Debug::ERROR) << " - " << filenames[fileIdx] << "\n";
            }
            Debug(Debug::ERROR) << "Please check your input files. Only files in fasta/fastq[.gz|bz2|zst] are supported\n";
            EXIT(EXIT_FAILURE);
        }
        size_t key = 0;
        size_t dataSize = 0;
        for (size_t i = 0; i < globalHistogram.size(); i++) {
            globalHistogram[i].startKey = key;
            key += globalHistogram[i].sequenceCount;
            dataSize += globalHistogram[i].sequenceCount * (globalHistogram[i].length + 2);
        }
        const size_t chunkSize = ceil(static_cast<double>(dataSize) / static_cast<double>(MMseqsNode::count));
        std::vector<size_t> dataFileSizes(2 * MMseqsNode::count, 0);
        size_t file = 0;
        size_t sumCharsAssigned = 0;
        std::string indexData;
        std::string headerIndexData;
        for (size_t i = 0; i < globalHistogram.size(); i++) {
            LengthBin &bin = globalHistogram[i];
            if (i > 0 && sumCharsAssigned >= (file + 1) * chunkSize) {
                file++;
            }
            sumCharsAssigned += bin.sequenceCount * (bin.length + 2);
            bin.file = file;
            bin.offset = dataFileSizes[file];
            bin.headerOffset = dataFileSizes[MMseqsNode::count + file];
            dataFileSizes[file] += bin.sequenceCount * (bin.length + 2);
            dataFileSizes[MMseqsNode::count + file] += bin.headerSize;
            indexData += SSTR(bin.startKey) + "\t" + SSTR(bin.offset) + "\t" + SSTR(bin.length + 2) + "\t" + SSTR(bin.sequenceCount) + "\t" + SSTR(bin.file) + "\n";
            headerIndexData += SSTR(bin.startKey) + "\t" + SSTR(bin.headerOffset) + "\t" + SSTR(bin.headerSize) + "\t" + SSTR(bin.sequenceCount) + "\t" + SSTR(bin.file) + "\n";
        }
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (size_t i = 0; i < dataFileSizes.size(); i++) {
            const std::string name = dataFile + (i < MMseqsNode::count ? "." : "_h.") + SSTR(i % MMseqsNode::count);
            FILE *handle = FileUtil::openAndDelete(name.c_str(), "w");
            if (ftruncate(fileno(handle), dataFileSizes[i]) != 0) {
                Debug(Debug::ERROR) << "Could not resize file " << name << "!\n";
                EXIT(EXIT_FAILURE);
            }
            if (fsync(fileno(handle)) != 0 || fclose(handle) != 0) {
                Debug(Debug::ERROR) << "Cannot close data file " << name << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        splitStarts(globalHistogram, nodeHistograms);
        for (unsigned int node = 0; node < MMseqsNode::count; node++) {
            std::string startsData;
            for (size_t i = 0; i < nodeHistograms[node].size(); i++) {
                const LengthBin &bin = nodeHistograms[node][i];
                startsData += SSTR(bin.length) + "\t" + SSTR(bin.sequenceCount) + "\t" + SSTR(bin.headerSize) + "\t" + SSTR(bin.startKey) + "\t" + SSTR(bin.offset) + "\t" + SSTR(bin.headerOffset) + "\t" + SSTR(bin.file) + "\n";
            }
            MMseqsNode::publishFile(startsPrefix + SSTR(node), startsData);
        }
        MMseqsNode::publishFile(dataFile + ".index", indexData);
        MMseqsNode::publishFile(dataFile + "_h.index", headerIndexData);
    }
    const std::string startsFile = startsPrefix + SSTR(MMseqsNode::rank);
    const std::vector<LengthBin> localNodeStarts = readLinfinityFile<LengthBin>(startsFile);
    FileUtil::remove(startsFile.c_str());
    for (size_t i = 0; i < pieceHistograms.size(); i++) {
        if (pieceHistograms[i].empty() == false && pieceHistograms[i].back().length == 0) {
            pieceHistograms[i].pop_back();
        }
    }
    splitStarts(localNodeStarts, pieceHistograms);

    std::vector<size_t> dataFileEnds(2 * MMseqsNode::count, 0);
    for (size_t i = 0; i < localNodeStarts.size(); i++) {
        const LengthBin &bin = localNodeStarts[i];
        dataFileEnds[bin.file] = bin.offset + bin.sequenceCount * (bin.length + 2);
        dataFileEnds[MMseqsNode::count + bin.file] = bin.headerOffset + bin.headerSize;
    }
    std::vector<int> dataFileFds(2 * MMseqsNode::count);
    for (size_t i = 0; i < dataFileFds.size(); i++) {
        const std::string name = dataFile + (i < MMseqsNode::count ? "." : "_h.") + SSTR(i % MMseqsNode::count);
        dataFileFds[i] = open(name.c_str(), O_WRONLY | O_CREAT, 0666);
        struct stat sb;
        if (dataFileFds[i] < 0 || fstat(dataFileFds[i], &sb) != 0 || static_cast<size_t>(sb.st_size) < dataFileEnds[i]) {
            Debug(Debug::ERROR) << "Cannot open data file " << name << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    const size_t activeThreads = std::max<size_t>(1, std::min<size_t>(par.threads, nodePieces.size()));
    size_t maxPieceSize = 0;
    for (size_t i = 0; i < pieceHistograms.size(); i++) {
        size_t pieceSize = 0;
        for (size_t j = 0; j < pieceHistograms[i].size(); j++) {
            pieceSize += pieceHistograms[i][j].sequenceCount * (pieceHistograms[i][j].length + 2) + pieceHistograms[i][j].headerSize;
        }
        maxPieceSize = std::max(maxPieceSize, pieceSize);
    }
    const size_t bufferSize = std::min(maxPieceSize, static_cast<size_t>(Util::computeMemory(0) * 0.95) / (2 * activeThreads));
#pragma omp parallel num_threads(activeThreads)
    {
        char *buffer = new(std::nothrow) char[bufferSize];
        Util::checkAllocation(buffer, "Can not allocate buffer memory in createdb");
        char *sortedBuffer = new(std::nothrow) char[bufferSize];
        Util::checkAllocation(sortedBuffer, "Can not allocate sortedBuffer memory in createdb");
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < nodePieces.size(); i++) {
            std::vector<LengthBin> lengthToBin = getLengthToBin(pieceHistograms[i]);
            size_t bufferPos = 0;
            KSeqWrapper *kseq = KSeqFactory(filenames[nodePieces[i].file].c_str(), nodePieces[i].start, nodePieces[i].end);
            while (kseq->ReadEntry()) {
                const KSeqWrapper::KSeqEntry &e = kseq->entry;
                if (e.sequence.l == 0) {
                    continue;
                }
                const size_t nameLength = strlen(e.name.s);
                const size_t entrySize = e.sequence.l + nameLength + 4;
                if (bufferPos + entrySize > bufferSize) {
                    if (entrySize > bufferSize) {
                        Debug(Debug::ERROR) << "Not enough memory for entry " << e.name.s << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    flushEntryBuffer(dataFileFds, dataFile, buffer, bufferPos, sortedBuffer, lengthToBin);
                    bufferPos = 0;
                }
                memcpy(buffer + bufferPos, e.sequence.s, e.sequence.l);
                std::replace(buffer + bufferPos, buffer + bufferPos + e.sequence.l, '\0', 'X');
                bufferPos += e.sequence.l;
                buffer[bufferPos++] = '\n';
                buffer[bufferPos++] = '\0';
                memcpy(buffer + bufferPos, e.name.s, nameLength);
                bufferPos += nameLength;
                buffer[bufferPos++] = '\n';
                buffer[bufferPos++] = '\0';
            }
            if (bufferPos > 0) {
                flushEntryBuffer(dataFileFds, dataFile, buffer, bufferPos, sortedBuffer, lengthToBin);
            }
            delete kseq;
        }
        delete[] buffer;
        delete[] sortedBuffer;
    }
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
    for (size_t i = 0; i < dataFileFds.size(); i++) {
        if (fsync(dataFileFds[i]) != 0 || close(dataFileFds[i]) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << dataFile << (i < MMseqsNode::count ? "." : "_h.") << (i % MMseqsNode::count) << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (par.compressHeaders) {
        MMseqsNode::publishFile(dataFile + ".written." + SSTR(MMseqsNode::rank), "");
        if (MMseqsNode::isMaster()) {
            MMseqsNode::gatherFiles(dataFile + ".written.");
            MMseqsNode::publishFile(dataFile + ".written", "");
        }
        MMseqsNode::waitForFile(dataFile + ".written");
        const std::string hdrDataFile = dataFile + "_h." + SSTR(MMseqsNode::rank);
        FILE *hdrFile = FileUtil::openFileOrDie(hdrDataFile.c_str(), "r", true);
        const size_t headerFileSize = FileUtil::getFileSize(hdrDataFile);
        const size_t frameSize = 1024 * 1024;
        const size_t frameBound = ZSTD_compressBound(frameSize);
        const size_t frameCount = (headerFileSize + frameSize - 1) / frameSize;
        const size_t batchSize = std::min(frameCount, static_cast<size_t>(4 * par.threads));
        char *frames = new(std::nothrow) char[batchSize * (frameSize + frameBound)];
        Util::checkAllocation(frames, "Can not allocate frame memory in createdb");
        std::vector<uint32_t> seekTable(2 * frameCount);
        const std::string tmpFile = hdrDataFile + ".tmp";
        FILE *zstdFile = FileUtil::openAndDelete(tmpFile.c_str(), "w");
        for (size_t batch = 0; batch < frameCount; batch += batchSize) {
            const size_t batchEnd = std::min(frameCount, batch + batchSize);
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (size_t frame = batch; frame < batchEnd; frame++) {
                char *input = frames + (frame - batch) * (frameSize + frameBound);
                const size_t size = std::min(frameSize, headerFileSize - frame * frameSize);
                for (size_t read = 0; read < size; ) {
                    const ssize_t count = pread(fileno(hdrFile), input + read, size - read, frame * frameSize + read);
                    if (count <= 0) {
                        Debug(Debug::ERROR) << "Can not read data file " << hdrDataFile << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    read += count;
                }
                const size_t compressedSize = ZSTD_compress(input + frameSize, frameBound, input, size, 3);
                if (ZSTD_isError(compressedSize)) {
                    Debug(Debug::ERROR) << "ZSTD_compress " << ZSTD_getErrorName(compressedSize) << "\n";
                    EXIT(EXIT_FAILURE);
                }
                seekTable[2 * frame] = compressedSize;
                seekTable[2 * frame + 1] = size;
            }
            for (size_t frame = batch; frame < batchEnd; frame++) {
                const char *output = frames + (frame - batch) * (frameSize + frameBound) + frameSize;
                if (fwrite(output, sizeof(char), seekTable[2 * frame], zstdFile) != seekTable[2 * frame]) {
                    Debug(Debug::ERROR) << "Can not write to data file " << tmpFile << "\n";
                    EXIT(EXIT_FAILURE);
                }
            }
        }
        delete[] frames;
        const uint32_t skippableHeader[2] = {ZSTD_MAGIC_SKIPPABLE_START | 0xE, static_cast<uint32_t>(seekTable.size() * sizeof(uint32_t) + 9)};
        const uint32_t footerFrames = frameCount;
        const char footerDescriptor = 0;
        const uint32_t footerMagic = 0x8F92EAB1;
        if (fwrite(skippableHeader, sizeof(uint32_t), 2, zstdFile) != 2
            || fwrite(seekTable.data(), sizeof(uint32_t), seekTable.size(), zstdFile) != seekTable.size()
            || fwrite(&footerFrames, sizeof(uint32_t), 1, zstdFile) != 1
            || fwrite(&footerDescriptor, sizeof(char), 1, zstdFile) != 1
            || fwrite(&footerMagic, sizeof(uint32_t), 1, zstdFile) != 1) {
            Debug(Debug::ERROR) << "Can not write to data file " << tmpFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fflush(zstdFile) != 0 || fsync(fileno(zstdFile)) != 0 || fclose(zstdFile) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << tmpFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(hdrFile) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << hdrDataFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (rename(tmpFile.c_str(), hdrDataFile.c_str()) != 0) {
            Debug(Debug::ERROR) << "Cannot rename " << tmpFile << " to " << hdrDataFile << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    MMseqsNode::publishFile(dataFile + ".done." + SSTR(MMseqsNode::rank), "");
    if (MMseqsNode::isMaster()) {
        MMseqsNode::gatherFiles(dataFile + ".done.");
        if (par.compressHeaders) {
            FileUtil::remove((dataFile + ".written").c_str());
        }
        std::string sourceData;
        for (size_t i = 0; i < filenames.size(); i++) {
            sourceData += SSTR(i) + "\t" + FileUtil::baseName(filenames[i]) + "\n";
        }
        MMseqsNode::publishFile(dataFile + ".source", sourceData);
        const int dbType = DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_AMINO_ACIDS, Parameters::DBTYPE_EXTENDED_LINFINITY);
        DBWriter::writeDbtypeFile(dataFile.c_str(), dbType, false);
        DBWriter::writeDbtypeFile((dataFile + "_h").c_str(), DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_GENERIC_DB, Parameters::DBTYPE_EXTENDED_LINFINITY), par.compressHeaders);
        Debug(Debug::INFO) << "Database type: " << Parameters::getDbTypeName(dbType) << "\n";
    }

    return EXIT_SUCCESS;
}

int createdb(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_VARIADIC, 0);
    if(par.gpu){
        par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_GPU;
        par.maskMode = 1;
    }else{
        par.maskMode = 0;
    }
    if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH) {
        MMseqsNode::init();
    }
    par.printParameters(command.cmd, argc, argv, *command.params);

    std::vector<std::string> filenames(par.filenames);
    std::string dataFile = filenames.back();
    filenames.pop_back();

    if (Util::endsWith(".tsv", filenames[0])) {
	    if (filenames.size() > 1) {
		    Debug(Debug::ERROR) << "Only one tsv file can be given\n";
		    EXIT(EXIT_FAILURE);
	    }
	    std::string tsv = filenames.back();
	    filenames.pop_back();

	    FILE* file = FileUtil::openFileOrDie(tsv.c_str(), "r", true);
	    char* line = NULL;
	    size_t len = 0;
	    ssize_t read;
	    while ((read = getline(&line, &len, file)) != -1) {
		    if (line[read - 1] == '\n') {
			    line[read - 1] = '\0';
			    read--;
		    }
		    filenames.push_back(line);
	    }
	    free(line);
	    fclose(file);
    }

    // consistent order
    SORT_SERIAL(filenames.begin(), filenames.end(), [](const std::string &a, const std::string &b) {
        return FileUtil::baseName(a) < FileUtil::baseName(b);
    });
    
    for (size_t i = 0; i < filenames.size(); i++) {
        if (FileUtil::directoryExists(filenames[i].c_str()) == true) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory\n";
            EXIT(EXIT_FAILURE);
        }
    }

    bool dbInput = false;
    if (FileUtil::fileExists(par.db1dbtype.c_str()) == true) {
        if (filenames.size() > 1) {
            Debug(Debug::ERROR) << "Only one database can be used with database input\n";
            EXIT(EXIT_FAILURE);
        }
        dbInput = true;
        par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
    }

    int dbType = -1;
    BaseMatrix * subMat = NULL;
    if (par.dbType == 2) {
        dbType = Parameters::DBTYPE_NUCLEOTIDES;
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 2.0, -0.0f);
    } else if(par.dbType == 1) {
        dbType = Parameters::DBTYPE_AMINO_ACIDS;
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.0f);
    }

    std::string indexFile = dataFile + ".index";
    if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && par.shuffleDatabase) {
        Debug(Debug::WARNING) << "Shuffle database cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_SOFT << "\n";
        Debug(Debug::WARNING) << "We recompute with --shuffle 0\n";
        par.shuffleDatabase = false;
    }

    if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_GPU && par.shuffleDatabase == false) {
        Debug(Debug::WARNING) << "Shuffle database cannot be turned off for --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_GPU << "\n";
        Debug(Debug::WARNING) << "We recompute with --shuffle 1\n";
        par.shuffleDatabase = true;
    }

    if ((par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT || par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH) && par.filenames[0] == "stdin") {
        Debug(Debug::WARNING) << "Stdin input cannot be combined with --createdb-mode " << par.createdbMode << "\n";
        Debug(Debug::WARNING) << "We recompute with --createdb-mode 1\n";
        par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
    }

    if ((par.maskMode == 1 || par.maskNrepeats > 0) && par.createdbMode != Parameters::SEQUENCE_SPLIT_MODE_GPU) {
        Debug(Debug::WARNING) << "Masking can only be used with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_GPU << "\n";
        Debug(Debug::WARNING) << "We recompute with --mask 0 or --mask-n-repeat 0\n";
        par.maskMode = false;
        par.maskNrepeats = 0;
    }

    const unsigned int shuffleSplits = par.shuffleDatabase ? 32 : 1;
    if ((par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT || par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH) && par.compressed) {
        Debug(Debug::WARNING) << "Compressed database cannot be combined with --createdb-mode " << par.createdbMode << "\n";
        Debug(Debug::WARNING) << "We recompute with --compressed 0\n";
        par.compressed = 0;
    }

    if (par.createdbMode != Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH && par.compressHeaders) {
        Debug(Debug::WARNING) << "Header compression can only be used with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << "\n";
        Debug(Debug::WARNING) << "We recompute with --compress-headers 0\n";
        par.compressHeaders = 0;
    }

    if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH && par.identifierOffset != 0) {
        Debug(Debug::WARNING) << "Identifier offset cannot be combined with --createdb-mode " << par.createdbMode << "\n";
        Debug(Debug::WARNING) << "We recompute with --id-offset 0\n";
        par.identifierOffset = 0;
    }

    if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH) {
        return createdbLinfinity(par, filenames, dataFile);
    }

    std::string hdrDataFile = dataFile + "_h";
    std::string hdrIndexFile = dataFile + "_h.index";

    DBKeyType entries_num = 0;
    const char newline = '\n';

    size_t sampleCount = 0;
    const size_t testForNucSequence = 100;
    size_t isNuclCnt = 0;
    Debug::Progress progress;
    std::vector<unsigned int>* sourceLookup = new std::vector<unsigned int>[shuffleSplits]();
    for (size_t i = 0; i < shuffleSplits; ++i) {
        sourceLookup[i].reserve(16384);
    }
    Debug(Debug::INFO) << "Converting sequences\n";

    std::string sourceFile = dataFile + ".source";

    redoComputation:
    FILE *source = fopen(sourceFile.c_str(), "w");
    if (source == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << sourceFile << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    DBWriter hdrWriter(hdrDataFile.c_str(), hdrIndexFile.c_str(), shuffleSplits, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    hdrWriter.open();
    DBWriter seqWriter(dataFile.c_str(), indexFile.c_str(), shuffleSplits, par.compressed, (dbType == -1) ? Parameters::DBTYPE_OMIT_FILE : dbType );
    seqWriter.open();
    size_t headerFileOffset = 0;
    size_t seqFileOffset = 0;

    size_t fileCount = filenames.size();
    DBReader<DBKeyType>* reader = NULL;
    if (dbInput == true) {
        reader = new DBReader<DBKeyType>(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_LOOKUP);
        reader->open(DBReader<DBKeyType>::LINEAR_ACCCESS);
        fileCount = reader->getSize();
    }

    // setup Sequence object pointer
    Sequence ** seqs = new Sequence * [par.threads];
    Masker ** masker = new Masker * [par.threads];
    masker[0] = NULL;
    const size_t BATCH_SIZE = par.threads * 10000;
    std::vector<std::pair<std::vector<char>, std::string>> batchEntries(BATCH_SIZE);
    size_t batchPos = 0;
    for (size_t fileIdx = 0; fileIdx < fileCount; fileIdx++) {
        size_t numEntriesInCurrFile = 0;
        std::string header;
        header.reserve(1024);

        std::string sourceName;
        if (dbInput == true) {
            DBKeyType dbKey = reader->getDbKey(fileIdx);
            size_t lookupId = reader->getLookupIdByKey(dbKey);
            sourceName = reader->getLookupEntryName(lookupId);
        } else {
            sourceName = FileUtil::baseName(filenames[fileIdx]);
        }
        char buffer[4096];
        size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\n", fileIdx, sourceName.c_str());
        size_t written = fwrite(buffer, sizeof(char), len, source);
        if (written != len) {
            Debug(Debug::ERROR) << "Cannot write to source file " << sourceFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        KSeqWrapper* kseq = NULL;
        if (dbInput == true) {
            kseq = new KSeqBuffer(reader->getData(fileIdx, 0), reader->getEntryLen(fileIdx) - 1);
        } else {
            kseq = KSeqFactory(filenames[fileIdx].c_str());
        }

        bool resetNotFile = par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && kseq->type != KSeqWrapper::KSEQ_FILE;
        if (resetNotFile) {
            Debug(Debug::WARNING) << "Only uncompressed fasta files can be used with --createdb-mode 0\n";
            Debug(Debug::WARNING) << "We recompute with --createdb-mode 1\n";
        }

        bool resetIncorrectNewline = false;
        if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT && kseq->type == KSeqWrapper::KSEQ_FILE) {
            // get last byte from filenames[fileIdx].c_str()
            FILE* fp = fopen(filenames[fileIdx].c_str(), "rb");
            if (fp == NULL) {
                Debug(Debug::ERROR) << "Cannot open file " << filenames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
            int res = fseek(fp, -1, SEEK_END);
            if (res != 0) {
                Debug(Debug::ERROR) << "Cannot seek at the end of file " << filenames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
            int lastChar = fgetc(fp);
            if (lastChar == EOF) {
                Debug(Debug::ERROR) << "Error reading from " << filenames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (fclose(fp) != 0) {
                Debug(Debug::ERROR) << "Error closing " << filenames[fileIdx] << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (lastChar != '\n') {
                Debug(Debug::WARNING) << "Last byte is not a newline. We recompute with --createdb-mode 1\n";
                resetIncorrectNewline = true;
            }
        }
        if (resetNotFile || resetIncorrectNewline) {
            par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
            progress.reset(SIZE_MAX);
            hdrWriter.close();
            seqWriter.close();
            delete kseq;
            if (fclose(source) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << sourceFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            for (size_t i = 0; i < shuffleSplits; ++i) {
                sourceLookup[i].clear();
            }
            goto redoComputation;
        }
        while (kseq->ReadEntry()) {
            progress.updateProgress();
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            if (e.name.l == 0) {
                Debug(Debug::ERROR) << "Fasta entry " << numEntriesInCurrFile << " is invalid\n";
                EXIT(EXIT_FAILURE);
            }

            DBKeyType id = static_cast<DBKeyType>(par.identifierOffset) + entries_num;
            if (dbType == -1) {
                // check for the first 10 sequences if they are nucleotide sequences
                if (sampleCount < 10 || (sampleCount % 100) == 0) {
                    if (sampleCount < testForNucSequence) {
                        size_t cnt = 0;
                        for (size_t i = 0; i < e.sequence.l; i++) {
                            switch (toupper(e.sequence.s[i])) {
                                case 'T':
                                case 'A':
                                case 'G':
                                case 'C':
                                case 'U':
                                case 'N':
                                    cnt++;
                                    break;
                            }
                        }
                        const float nuclDNAFraction = static_cast<float>(cnt) / static_cast<float>(e.sequence.l);
                        if (nuclDNAFraction > 0.9) {
                            isNuclCnt += true;
                        }
                    }
                    sampleCount++;
                }
            }

            if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT) {
                if (e.newlineCount != 1) {
                    if (e.newlineCount == 0) {
                        Debug(Debug::WARNING) << "Fasta entry " << numEntriesInCurrFile << " has no newline character\n";
                    } else if (e.newlineCount > 1) {
                        Debug(Debug::WARNING) << "Multiline fasta can not be combined with --createdb-mode 1\n";
                    }
                    Debug(Debug::WARNING) << "We recompute with --createdb-mode 0\n";
                    par.createdbMode = Parameters::SEQUENCE_SPLIT_MODE_HARD;
                    progress.reset(SIZE_MAX);
                    hdrWriter.close();
                    seqWriter.close();
                    delete kseq;
                    if (fclose(source) != 0) {
                        Debug(Debug::ERROR) << "Cannot close file " << sourceFile << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    for (size_t i = 0; i < shuffleSplits; ++i) {
                        sourceLookup[i].clear();
                    }
                    goto redoComputation;
                }
            }

            header.append(e.name.s, e.name.l);
            if (e.comment.l > 0) {
                header.append(" ", 1);
                header.append(e.comment.s, e.comment.l);
            }

            std::string headerId = Util::parseFastaHeader(header.c_str());
            if (headerId.empty()) {
                // An identifier is necessary for these two cases, so we should just give up
                Debug(Debug::WARNING) << "Cannot extract identifier from entry " << numEntriesInCurrFile << "\n";
            }
            header.push_back('\n');

            // Finally write down the entry
            unsigned int splitIdx = id % shuffleSplits;
            sourceLookup[splitIdx].emplace_back(fileIdx);
            if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT) {
                // +2 to emulate the \n\0
                hdrWriter.writeIndexEntry(id, headerFileOffset + e.headerOffset, (e.sequenceOffset-e.headerOffset)+1, 0);
                seqWriter.writeIndexEntry(id, seqFileOffset + e.sequenceOffset, e.sequence.l+2, 0);
            } else if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_HARD) {
                hdrWriter.writeData(header.c_str(), header.length(), id, splitIdx);
                seqWriter.writeStart(splitIdx);
                seqWriter.writeAdd(e.sequence.s, e.sequence.l, splitIdx);
                seqWriter.writeAdd(&newline, 1, splitIdx);
                seqWriter.writeEnd(id, splitIdx, true);
            } else {
                // Add to batch, allocated sequence and masker and genrate submat on the fly
                batchEntries[batchPos].first.clear();
                //batchEntries[batchPos].first.reserve(e.sequence.l);
                for(size_t j = 0; j < e.sequence.l; j++){
                    batchEntries[batchPos].first.push_back(e.sequence.s[j]);
                }
//                batchEntries[batchPos].first.push_back('\n');
                //std::copy_n(e.sequence.s, e.sequence.l, batchEntries[batchPos].first.begin());
                batchEntries[batchPos].second.clear();
                batchEntries[batchPos].second.append(header.begin(), header.end());
                batchPos++;
                if(batchPos == BATCH_SIZE){
                    if(subMat == NULL){
                        if (isNuclCnt == sampleCount) {
                            subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 2.0, -0.0f);
                            dbType = Parameters::DBTYPE_NUCLEOTIDES;
                        } else{
                            subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.0f);
                            dbType = Parameters::DBTYPE_AMINO_ACIDS;
                        }
                    }
                    processSeqBatch(par, seqWriter, hdrWriter, subMat, dbType, masker, seqs,
                                    id - (batchPos - 1), batchEntries, batchPos, shuffleSplits);
                    batchPos = 0;
                }
            }

            entries_num++;
            numEntriesInCurrFile++;
            header.clear();
        }

        if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_GPU && batchPos > 0){
            if(subMat == NULL){
                if (isNuclCnt == sampleCount) {
                    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 2.0, -0.0f);
                    dbType = Parameters::DBTYPE_NUCLEOTIDES;
                } else {
                    subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.0f);
                    dbType = Parameters::DBTYPE_AMINO_ACIDS;
                }
            }
            processSeqBatch(par, seqWriter, hdrWriter, subMat, dbType, masker, seqs,
                            (par.identifierOffset + entries_num) - (batchPos - 1), batchEntries, batchPos, shuffleSplits);
        }

        if (numEntriesInCurrFile == 0) {
            Debug(Debug::WARNING) << "File " << sourceName << " is empty or invalid and was ignored\n";
        }

        delete kseq;
        if (filenames.size() > 1 && par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT) {
            size_t fileSize = FileUtil::getFileSize(filenames[fileIdx].c_str());
            headerFileOffset += fileSize;
            seqFileOffset += fileSize;
        }
    }
    Debug(Debug::INFO) << "\n";
    if (fclose(source) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << sourceFile << "\n";
        EXIT(EXIT_FAILURE);
    }

    // sort
    bool gpuCompatibleDB = false;
    Timer timer;
    if(par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_GPU){
        gpuCompatibleDB = true;
        hdrWriter.closeFiles();
        seqWriter.closeFiles();
        for(unsigned int i = 0; i < shuffleSplits; i++){
            sortWithIndex(seqWriter.getDataFileNames()[i], seqWriter.getIndexFileNames()[i],
                          hdrWriter.getDataFileNames()[i], hdrWriter.getIndexFileNames()[i]);
        }
        Debug(Debug::INFO) << "Sort single files in " << timer.lap() << "\n";
        std::string lookupFile = dataFile + ".lookup";

        timer.reset();
        mergeSequentialByJointIndex(seqWriter.getDataFileNames(), seqWriter.getIndexFileNames(),
                                    hdrWriter.getDataFileNames(), hdrWriter.getIndexFileNames(),
                                    seqWriter.getDataFileName(), seqWriter.getIndexFileName(),
                                    hdrWriter.getDataFileName(), hdrWriter.getIndexFileName(),
                                    lookupFile.c_str(), sourceLookup, entries_num, shuffleSplits);
        Debug(Debug::INFO) << "Merge all files " << timer.lap() << "\n";
        hdrWriter.clearMemory();
        seqWriter.clearMemory();
        //TODO we cannot have compressed seq. dbs anymore if wanted ot have GPU compatible dbs
    } else {
        hdrWriter.close(true, false);
        seqWriter.close(true, false);
        if (par.shuffleDatabase == true) {
            DBWriter::createRenumberedDB(dataFile, indexFile, "", "", DBReader<DBKeyType>::LINEAR_ACCCESS);
            DBWriter::createRenumberedDB(hdrDataFile, hdrIndexFile, "", "", DBReader<DBKeyType>::LINEAR_ACCCESS);
        }
        if (par.createdbMode == Parameters::SEQUENCE_SPLIT_MODE_SOFT) {
            if (filenames.size() == 1) {
                FileUtil::symlinkAbs(filenames[0], dataFile);
                FileUtil::symlinkAbs(filenames[0], hdrDataFile);
            } else {
                for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
                    FileUtil::symlinkAbs(filenames[fileIdx], dataFile + "." + SSTR(fileIdx));
                    FileUtil::symlinkAbs(filenames[fileIdx], hdrDataFile + "." + SSTR(fileIdx));
                }
            }
        }
        if (par.writeLookup == true) {
            DBReader<DBKeyType> readerHeader(hdrDataFile.c_str(), hdrIndexFile.c_str(), 1, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX);
            readerHeader.open(DBReader<DBKeyType>::NOSORT);
            // create lookup file
            std::string lookupFile = dataFile + ".lookup";
            FILE* file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
            std::string buffer;
            buffer.reserve(2048);
            size_t splitIdx = 0;
            size_t splitCounter = 0;
            DBReader<DBKeyType>::LookupEntry entry;
            for (size_t id = 0; id < readerHeader.getSize(); id++) {
                size_t splitSize = sourceLookup[splitIdx].size();
                if (splitSize == 0 || splitCounter > sourceLookup[splitIdx].size() - 1) {
                    splitIdx++;
                    splitCounter = 0;
                }
                char *header = readerHeader.getData(id, 0);
                entry.id = id;
                entry.entryName = Util::parseFastaHeader(header);
                if (entry.entryName.empty()) {
                    Debug(Debug::WARNING) << "Cannot extract identifier from entry " << entries_num << "\n";
                }
                entry.fileNumber = sourceLookup[splitIdx][splitCounter];
                readerHeader.lookupEntryToBuffer(buffer, entry);
                int written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
                if (written != (int)buffer.size()) {
                    Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                    EXIT(EXIT_FAILURE);
                }
                buffer.clear();
                splitCounter++;
            }
            if (fclose(file) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            readerHeader.close();
        }
    }

    if (dbType == -1) {
        if (isNuclCnt == sampleCount) {
            dbType = Parameters::DBTYPE_NUCLEOTIDES;
        } else {
            dbType = Parameters::DBTYPE_AMINO_ACIDS;
        }
    }
    if(gpuCompatibleDB){
        dbType = DBReader<DBKeyType>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_GPU);
    }
    DBWriter::writeDbtypeFile(seqWriter.getDataFileName(), dbType ,par.compressed);
    DBWriter::writeDbtypeFile(hdrWriter.getDataFileName(), Parameters::DBTYPE_GENERIC_DB, par.compressed);

    Debug(Debug::INFO) << "Database type: " << Parameters::getDbTypeName(dbType) << "\n";
    if (dbInput == true) {
        reader->close();
        delete reader;
    }

    if (entries_num == 0) {
        Debug(Debug::ERROR) << "The input files have no entry: ";
        for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
            Debug(Debug::ERROR) << " - " << filenames[fileIdx] << "\n";
        }
        Debug(Debug::ERROR) << "Please check your input files. Only files in fasta/fastq[.gz|bz2|zst] are supported\n";
        EXIT(EXIT_FAILURE);
    }

    if(masker[0] != NULL){
        for(int i = 0; i < par.threads; i++){
            delete masker[i];
            delete seqs[i];
        }
    }
    delete [] masker;
    delete [] seqs;

    if (subMat != NULL) {
        delete subMat;
    }

    delete[] sourceLookup;

    return EXIT_SUCCESS;
}
