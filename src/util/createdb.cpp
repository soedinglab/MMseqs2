/*
 * createdb
 * written by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.
 * modified by Maria Hauser <mhauser@genzentrum.lmu.de> (splitting into sequences/headers databases)
 * modified by Milot Mirdita <milot@mirdita.de>
 */

#include <cstdio>

#include <map>
#include <fstream>
#include <unistd.h>
#include <math.h>
#include <itoa.h>
#include <random>

#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"
#include "KSeqWrapper.h"

int createdb(int argn, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argn, argv, command, 2, true, Parameters::PARSE_VARIADIC);

    if (par.maxSeqLen == Parameters::MAX_SEQ_LEN) {
        par.maxSeqLen = Parameters::MAX_SEQ_LEN - 1;
    }

    std::vector<std::string> filenames(par.filenames);

    std::string data_filename = filenames.back();
    filenames.pop_back();
    std::string index_filename = data_filename + ".index";

    std::string data_filename_hdr(data_filename);
    data_filename_hdr.append("_h");

    std::string index_filename_hdr(data_filename);
    index_filename_hdr.append("_h.index");

    std::string lookupFileName(data_filename);
    lookupFileName.append(".lookup");
    FILE* lookupFile = fopen(lookupFileName.c_str(), "w");
    if(lookupFile == NULL) {
        Debug(Debug::ERROR) << "Could not open " << lookupFileName << " for writing.";
        EXIT(EXIT_FAILURE);
    }

    for(size_t i = 0; i < filenames.size(); i++){
        if(FileUtil::fileExists(filenames[i].c_str())==false){
            Debug(Debug::ERROR) << "File " << filenames[i] << " does not exist.\n";
            EXIT(EXIT_FAILURE);
        }
        if(FileUtil::directoryExists(filenames[i].c_str())==true){
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    DBWriter out_writer(data_filename.c_str(), index_filename.c_str());
    DBWriter out_hdr_writer(data_filename_hdr.c_str(), index_filename_hdr.c_str());
    out_writer.open();
    out_hdr_writer.open();

    DBWriter *set_writer = NULL;
    if (par.clusterDB) {
        std::string cluster_filename(data_filename);
        cluster_filename.append("_member_lookup");
        std::string cluster_index_filename(data_filename);
        cluster_index_filename.append("_member_lookup.index");
        set_writer = new DBWriter(cluster_filename.c_str(), cluster_index_filename.c_str());
        set_writer->open();
    }

    unsigned int entries_num = 1;
    size_t count = 0;
    size_t sampleCount = 0;

    const size_t testForNucSequence = 100;
    size_t isNuclCnt = 0;

    for (size_t i = 0; i < filenames.size(); i++) {
        std::string splitHeader;
        splitHeader.reserve(1024);
        std::string header;
        header.reserve(1024);
        std::string splitId;
        splitId.reserve(1024);
        KSeqWrapper *kseq = KSeqFactory(filenames[i].c_str());
        std::string setBuffer;
        char keyBuffer[255];
        while (kseq->ReadEntry()) {
            Debug::printProgress(count);
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            if (e.name.length() == 0) {
                Debug(Debug::ERROR) << "Fasta entry: " << entries_num << " is invalid.\n";
                EXIT(EXIT_FAILURE);
            }

            size_t splitCnt = 1;
            if (par.splitSeqByLen == true) {
                splitCnt = (size_t) ceilf(static_cast<float>(e.sequence.length()) / static_cast<float>(par.maxSeqLen));
            }

            // header
            header.append(e.name);
            if (e.comment.length() > 0) {
                header.append(" ", 1);
                header.append(e.comment);
            }

            std::string headerId = Util::parseFastaHeader(header);
            if (headerId == "") {
                // An identifier is necessary for these two cases, so we should just give up
                Debug(Debug::WARNING) << "Could not extract identifier from entry " << entries_num << ".\n";

            }
            for (size_t split = 0; split < splitCnt; split++) {
                splitId.append(headerId);
                if (splitCnt > 1) {
                    splitId.append("_");
                    splitId.append(SSTR(split));
                }

                unsigned int id = par.identifierOffset + entries_num;

                // For split entries replace the found identifier by identifier_splitNumber
                // Also add another hint that it was split to the end of the header
                splitHeader.append(header);
                if (par.splitSeqByLen == true && splitCnt > 1) {
                    if (headerId != "") {
                        size_t pos = splitHeader.find(headerId);
                        if (pos != std::string::npos) {
                            splitHeader.erase(pos, headerId.length());
                            splitHeader.insert(pos, splitId);
                        }
                    }
                    splitHeader.append(" Split=");
                    splitHeader.append(SSTR(split));
                }

                // space is needed for later parsing
                splitHeader.append(" ", 1);
                splitHeader.append("\n");

                // Finally write down the entry
                out_hdr_writer.writeData(splitHeader.c_str(), splitHeader.length(), id);
                splitHeader.clear();
                splitId.clear();

                // sequence
                const std::string &sequence = e.sequence;
                // check for the first 10 sequences if they are nucleotide sequences
                if((count % 100) == 0){
                    if(sampleCount < testForNucSequence){
                        size_t cnt=0;
                        for(size_t i = 0; i < sequence.size(); i++){
                            switch(toupper(sequence[i]))
                            {
                                case 'T':
                                case 'A':
                                case 'G':
                                case 'C':
                                case 'N': cnt++;
                                break;
                            }
                        }
                        if(cnt == sequence.size()){
                            isNuclCnt += true;
                        }
                    }
                    sampleCount++;
                }

                if (par.splitSeqByLen) {
                    size_t len = std::min(par.maxSeqLen, sequence.length() - split * par.maxSeqLen);
                    std::string splitString(sequence.c_str() + split * par.maxSeqLen, len);
                    splitString.append("\n");
                    out_writer.writeData(splitString.c_str(), splitString.length(), id);
                } else {
                    std::string seqWithLineRet(sequence);
                    seqWithLineRet.append("\n");
                    out_writer.writeData(seqWithLineRet.c_str(), seqWithLineRet.length(), id);
                }

                if (set_writer != NULL) {
                    Itoa::u32toa_sse2(id, keyBuffer);
                    setBuffer.append(keyBuffer);
                    setBuffer.append("\n");
                }
                entries_num++;
                count++;
            }
            header.clear();
        }
        if (set_writer != NULL) {
            set_writer->writeData(setBuffer.c_str(), setBuffer.size(), i);
        }
        setBuffer.clear();
        delete kseq;
    }

    int dbType = Sequence::AMINO_ACIDS;
    if (isNuclCnt == sampleCount || isNuclCnt == testForNucSequence) {
        dbType = Sequence::NUCLEOTIDES;
    }
    out_hdr_writer.close();
    out_writer.close(dbType);
    if(par.clusterDB){
        set_writer->close();
    }


    DBReader<unsigned int> readerHeader(out_hdr_writer.getDataFileName(), out_hdr_writer.getIndexFileName());
    readerHeader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::Index *indexHeader = readerHeader.getIndex();
    unsigned int *lengthHeader = readerHeader.getSeqLens();

    if(!par.clusterDB) {
        // shuffle data
        DBReader<unsigned int> readerSequence(out_writer.getDataFileName(), out_writer.getIndexFileName());
        readerSequence.open(DBReader<unsigned int>::NOSORT);
        readerSequence.readMmapedDataInMemory();
        DBReader<unsigned int>::Index *indexSequence = readerSequence.getIndex();
        unsigned int *lengthSequence = readerSequence.getSeqLens();

        unsigned int *perm = new unsigned int[readerSequence.getSize()];
        std::default_random_engine generator(0);
        std::uniform_int_distribution<unsigned int> distribution(0, readerSequence.getSize() - 1);
        for (unsigned int n = 0; n < readerSequence.getSize(); n++) {
            perm[n] = n;
        }
        for (unsigned int n = 0; n < readerSequence.getSize(); n++) {
            unsigned int n_new = distribution(generator);
            std::swap(indexSequence[n_new], indexSequence[n]);
            std::swap(lengthSequence[n_new], lengthSequence[n]);
            std::swap(indexHeader[n_new], indexHeader[n]);
            std::swap(lengthHeader[n_new], lengthHeader[n]);
        }
        delete[] perm;
        DBWriter out_writer_shuffeled(data_filename.c_str(), index_filename.c_str());
        out_writer_shuffeled.open();
        for (unsigned int n = 0; n < readerSequence.getSize(); n++) {
            unsigned int id = par.identifierOffset + n;
            const char *data = readerSequence.getData() + indexSequence[n].offset;
            out_writer_shuffeled.writeData(data, lengthSequence[n] - 1, id);
        }
        readerSequence.close();
        out_writer_shuffeled.close(dbType);

        DBWriter out_hdr_writer_shuffeled(data_filename_hdr.c_str(), index_filename_hdr.c_str());
        out_hdr_writer_shuffeled.open();
        readerHeader.readMmapedDataInMemory();
        char lookupBuffer[32768];
        for (unsigned int n = 0; n < readerHeader.getSize(); n++) {
            unsigned int id = par.identifierOffset + n;
            const char *data = readerHeader.getData() + indexHeader[n].offset;
            std::string splitId = Util::parseFastaHeader(data);
            char *tmpBuff = Itoa::u32toa_sse2(id, lookupBuffer);
            *(tmpBuff - 1) = '\t';
            fwrite(lookupBuffer, sizeof(char), tmpBuff - lookupBuffer, lookupFile);
            fwrite(splitId.c_str(), sizeof(char), splitId.length(), lookupFile);
            char newline = '\n';
            fwrite(&newline, sizeof(char), 1, lookupFile);
            out_hdr_writer_shuffeled.writeData(data, lengthHeader[n] - 1, id);
        }
        out_hdr_writer_shuffeled.close();

    } else {
        readerHeader.readMmapedDataInMemory();
        char lookupBuffer[32768];
        for (unsigned int n = 0; n < readerHeader.getSize(); n++) {
            unsigned int id = par.identifierOffset + n + 1;
            const char *data = readerHeader.getData() + indexHeader[n].offset;
            std::string splitId = Util::parseFastaHeader(data);
            char *tmpBuff = Itoa::u32toa_sse2(id, lookupBuffer);
            *(tmpBuff - 1) = '\t';
            fwrite(lookupBuffer, sizeof(char), tmpBuff - lookupBuffer, lookupFile);
            fwrite(splitId.c_str(), sizeof(char), splitId.length(), lookupFile);
            char newline = '\n';
            fwrite(&newline, sizeof(char), 1, lookupFile);
        }
    }
    readerHeader.close();
    fclose(lookupFile);

    return EXIT_SUCCESS;
}




