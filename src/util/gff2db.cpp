#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MemoryMapped.h"
#include "Orf.h"
#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

int gff2db(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string outDb = par.filenames.back();
    par.filenames.pop_back();
    std::string seqDb = par.filenames.back();
    par.filenames.pop_back();

    DBReader<unsigned int> reader(seqDb.c_str(), (seqDb + ".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_LOOKUP_REV);
    reader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> headerReader((seqDb + "_h").c_str(), (seqDb + "_h.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::NOSORT);

    std::string outDbIndex = outDb + ".index";
    DBWriter writer(outDb.c_str(), outDbIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    writer.open();
    std::string outHdr = outDb + "_h";
    std::string outHdrIndex = outDb + "_h.index";
    DBWriter headerWriter(outHdr.c_str(), outHdrIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();
    std::string outLookup = outDb + ".lookup";
    std::string outLookupIndex = outDb + ".lookup.index";
    DBWriter lookupWriter(outLookup.c_str(), outLookupIndex.c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    lookupWriter.open();

    FILE *source = FileUtil::openAndDelete((outDb + ".source").c_str(), "w");
    for (size_t i = 0; i < par.filenames.size(); ++i) {
        if (fprintf(source, "%zu\t%s\n", i, FileUtil::baseName(par.filenames[i]).c_str()) < 0) {
            Debug(Debug::ERROR) << "Cannot write to file " << outDb << ".source\n";
            return EXIT_FAILURE;
        }
    }
    if (fclose(source) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << outDb << ".source\n";
        return EXIT_FAILURE;
    }

    const std::vector<std::string> features = Util::split(par.gffType, ",");
    if (features.empty()) {
        Debug(Debug::WARNING) << "No feature types given. All features will be extracted\n";
    }
    std::vector<size_t> featureCount(features.size(), 0);

    if (par.filenames.size() < reader.getSize()) {
        Debug(Debug::WARNING) << "Not enough GFF files are provided. Some results might be omitted\n";
    }

    unsigned int entries_num = 0;
    Debug::Progress progress(par.filenames.size());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        char buffer[32768];
        const char* fields[255];

        std::string revStr;
        revStr.reserve(par.maxSeqLen + 1);

        std::vector<size_t> localFeatureCount(features.size(), 0);

#pragma omp for schedule(dynamic, 1) nowait
        for (size_t i = 0; i < par.filenames.size(); ++i) {
            progress.updateProgress();
            MemoryMapped file(par.filenames[i], MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
            if (!file.isValid()) {
                Debug(Debug::ERROR) << "Could not open GFF file " << par.filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
            char *data = (char *) file.getData();
            char* end = data + file.mappedSize();
            size_t idx = 0;
            while (data < end && *data != '\0') {
                // line is a comment or empty
                if (*data == '#' || *data == '\n') {
                    data = Util::skipLine(data);
                    continue;
                }

                const size_t columns = Util::getFieldsOfLine(data, fields, 255);
                data = Util::skipLine(data);
                if (columns < 9) {
                    Debug(Debug::WARNING) << "Not enough columns in GFF file\n";
                    continue;
                }

                if (features.empty() == false) {
                    bool shouldSkip = true;
                    std::string type(fields[2], fields[3] - fields[2] - 1);
                    for (size_t i = 0; i < features.size(); ++i) {
                        if (type.compare(features[i]) == 0) {
                            localFeatureCount[i]++;
                            shouldSkip = false;
                            break;
                        }
                    }
                    if (shouldSkip) {
                        continue;
                    }
                }
                size_t start = Util::fast_atoi<size_t>(fields[3]);
                size_t end = Util::fast_atoi<size_t>(fields[4]);
                if (start == end) {
                    Debug(Debug::WARNING) << "Invalid sequence length in line " << idx << "\n";
                    continue;
                }
                std::string strand(fields[6], fields[7] - fields[6] - 1);
                std::string name(fields[0], fields[1] - fields[0] - 1);
                size_t lookupId = reader.getLookupIdByAccession(name);
                if (lookupId == SIZE_MAX) {
                    Debug(Debug::ERROR) << "GFF entry not found in database lookup: " << name << "\n";
                    EXIT(EXIT_FAILURE);
                }
                unsigned int lookupKey = reader.getLookupKey(lookupId);
                size_t seqId = reader.getId(lookupKey);
                if (seqId == UINT_MAX) {
                    Debug(Debug::ERROR) << "GFF entry not found in sequence database: " << name << "\n";
                    EXIT(EXIT_FAILURE);
                }

                unsigned int key = __sync_fetch_and_add(&(entries_num), 1);
                size_t bufferLen;
                if (strand == "+") {
                    bufferLen = Orf::writeOrfHeader(buffer, lookupKey, start, end, 0, 0);
                } else {
                    bufferLen = Orf::writeOrfHeader(buffer, lookupKey, end, start, 0, 0);
                }
                headerWriter.writeData(buffer, bufferLen, key, thread_idx);

                char *seq = reader.getData(seqId, thread_idx);
                //check the strand instead of end - start 
                ssize_t length = end - start + 1;
                writer.writeStart(thread_idx);
                if (strand == "+") {
                    size_t len = snprintf(buffer, sizeof(buffer), "%u\t%s_%zu_%zu_%zu\t%zu\n", key, name.c_str(), idx, start, end, i);
                    lookupWriter.writeData(buffer, len, thread_idx, false, false);
                    writer.writeAdd(seq + start - 1 , length, thread_idx);
                } else {
                    size_t len = snprintf(buffer, sizeof(buffer), "%u\t%s_%zu_%zu_%zu\t%zu\n", key, name.c_str(), idx, end, start, i);
                    lookupWriter.writeData(buffer, len, thread_idx, false, false);
                    for (size_t i = end - 1 ; i >= end - length; i--) {
                        revStr.append(1, Orf::complement(seq[i]));
                    }
                    writer.writeAdd(revStr.c_str(), revStr.size(), thread_idx);
                    revStr.clear();
                }
                writer.writeAdd("\n", 1, thread_idx);
                writer.writeEnd(key, thread_idx);
                idx++;
            }
            file.close();
        }
#pragma omp critical
        for (size_t i = 0; i < features.size(); ++i) {
            featureCount[i] += localFeatureCount[i];
        }
    }
    lookupWriter.close(true);
    FileUtil::remove(lookupWriter.getIndexFileName());
    headerWriter.close(true);
    writer.close(true);
    headerReader.close();
    reader.close();
    if (Debug::debugLevel >= Debug::INFO && features.size() > 0) {
        Debug(Debug::INFO) << "Found these feature types and counts:\n";
        for (size_t i = 0; i < features.size(); ++i) {
            Debug(Debug::INFO) << " - " << features[i] << ": " << featureCount[i] << "\n";
        }
    } else {
        Debug(Debug::INFO) << (entries_num + 1) << " features were extracted\n";
    }

    if (par.filenames.size() > 1 && par.threads > 1) {
    // make identifiers stable
#pragma omp parallel
        {
#pragma omp single
            {
#pragma omp task
                {
                    DBWriter::createRenumberedDB(outHdr, outHdrIndex, "", "");
                }

#pragma omp task
                {
                    DBWriter::createRenumberedDB(outDb, outDbIndex, outDb, outDbIndex);
                }
            }
        }
    }

    return EXIT_SUCCESS;
}



