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

    DBWriter writer(outDb.c_str(), (outDb + ".index").c_str(), par.threads, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    writer.open();
    DBWriter headerWriter((outDb + "_h").c_str(), (outDb + "_h.index").c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();
    DBWriter lookupWriter((outDb + ".lookup").c_str(), (outDb + ".lookup.index").c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    lookupWriter.open();

    FILE *source = FileUtil::openAndDelete((outDb + ".source").c_str(), "w");
    for (size_t i = 0; i < par.filenames.size(); ++i) {
        if (fprintf(source, "%zu\t%s\n", i, FileUtil::baseName(par.filenames[i]).c_str()) < 0) {
            Debug(Debug::ERROR) << "Cannot write to file " << outDb << ".source\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (fclose(source) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << outDb << ".source\n";
        EXIT(EXIT_FAILURE);
    }

    bool shouldCompareType = par.gffType.length() > 0;
    if (par.filenames.size() < reader.getSize()) {
        Debug(Debug::WARNING) << "Not enough GFF files are provided. Some results might be omitted. \n";
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

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < par.filenames.size(); ++i) {
            progress.updateProgress();
            MemoryMapped file(par.filenames[i], MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
            if (!file.isValid()) {
                Debug(Debug::ERROR) << "Could not open GFF file " << par.filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
            char *data = (char *) file.getData();
            size_t idx = 0;
            while (*data != '\0') {
                // line is a comment
                if (*data == '#') {
                    data = Util::skipLine(data);
                    continue;
                }
                //TODO: multi-word column 2 could exist
                const size_t columns = Util::getWordsOfLine(data, fields, 255);
                data = Util::skipLine(data);
                if (columns < 9) {
                    Debug(Debug::WARNING) << "Not enough columns in GFF file\n";
                    continue;
                }

                if (shouldCompareType) {
                    std::string type(fields[2], fields[3] - fields[2] - 1);
                    if (type.compare(par.gffType) != 0) {
                        continue;
                    }
                }

                unsigned int key = __sync_fetch_and_add(&(entries_num), 1);
                size_t start = Util::fast_atoi<size_t>(fields[3]);
                size_t end = Util::fast_atoi<size_t>(fields[4]);
                if (start == end) {
                    Debug(Debug::WARNING) << "Invalid sequence length in line " << key << "!\n";
                    continue;
                }
                std::string strand(fields[6], fields[7] - fields[6] - 1);
                std::string name(fields[0], fields[1] - fields[0] - 1);
                size_t lookupId = reader.getLookupIdByAccession(name);
                if (lookupId == SIZE_MAX) {
                    Debug(Debug::ERROR) << "GFF entry not found in database lookup: " << name << "!\n";
                    EXIT(EXIT_FAILURE);
                }
                unsigned int lookupKey = reader.getLookupKey(lookupId);

                idx++;
                size_t bufferLen;
                if (strand == "+") {
                    bufferLen = Orf::writeOrfHeader(buffer, lookupKey, start, end, 0, 0);
                } else {
                    bufferLen = Orf::writeOrfHeader(buffer, lookupKey, end, start, 0, 0);
                }
                headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                

                size_t seqId = reader.getId(lookupKey);
                if (seqId == UINT_MAX) {
                    Debug(Debug::ERROR) << "GFF entry not found in sequence database: " << name << "!\n";
                    EXIT(EXIT_FAILURE);
                }

                //check the strand instead of end - start 
                ssize_t length = end - start + 1;
                char *seq = reader.getData(seqId, thread_idx);

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
                    if (revStr.size() % 3 != 0) {
                        Debug(Debug::WARNING) << "The nucleotide sequence length is not divisible by 3!\n";
                    }
                    writer.writeAdd(revStr.c_str(), revStr.size(), thread_idx);
                    revStr.clear();
                }
                writer.writeAdd("\n", 1, thread_idx);
                writer.writeEnd(key, thread_idx);
            }
            file.close();
        }
    }
    lookupWriter.close(true);
    FileUtil::remove(lookupWriter.getIndexFileName());
    headerWriter.close(true);
    writer.close(true);
    headerReader.close();
    reader.close();

    return EXIT_SUCCESS;
}



