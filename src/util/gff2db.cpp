#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MemoryMapped.h"
#include "Orf.h"

int gff2db(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    MemoryMapped file(par.db1, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (!file.isValid()) {
        Debug(Debug::ERROR) << "Could not open GFF file " << par.db1 << "\n";
        EXIT(EXIT_FAILURE);
    }
    char *data = (char *) file.getData();

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_LOOKUP_REV);
    reader.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> headerReader(par.hdr2.c_str(), par.hdr2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    headerReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    writer.open();
    DBWriter headerWriter(par.hdr3.c_str(), par.hdr3Index.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    bool shouldCompareType = par.gffType.length() > 0;

    Debug::Progress progress;
    unsigned int entries_num = 0;
    char buffer[1024];
    const char* fields[255];

    std::string revStr;
    revStr.reserve(par.maxSeqLen + 1);

    while (*data != '\0') {
        progress.updateProgress();

        // line is a comment
        if (*data == '#') {
            data = Util::skipLine(data);
            continue;
        }

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

        size_t start = Util::fast_atoi<size_t>(fields[3]);
        size_t end = Util::fast_atoi<size_t>(fields[4]);
        if (start == end) {
            Debug(Debug::WARNING) << "Invalid sequence length in line " << entries_num << "!\n";
            continue;
        }

        std::string name(fields[0], fields[1] - fields[0] - 1);
        size_t lookupId = reader.getLookupIdByAccession(name);
        if (lookupId == SIZE_MAX) {
            Debug(Debug::ERROR) << "GFF entry not found in database lookup: " << name << "!\n";
            return EXIT_FAILURE;
        }
        unsigned int key = reader.getLookupKey(lookupId);

        size_t headerId = headerReader.getId(key);
        if (headerId == UINT_MAX) {
            Debug(Debug::ERROR) << "GFF entry not found in header database: " << name << "!\n";
            return EXIT_FAILURE;
        }
        unsigned int id = par.identifierOffset + entries_num;

        headerWriter.writeStart(0);
        headerWriter.writeAdd(headerReader.getData(headerId, 0), std::max(headerReader.getSeqLen(headerId), (size_t)2) - 2, 0);
        int len = snprintf(buffer, 1024, " %zu-%zu\n", start, end);
        headerWriter.writeAdd(buffer, len, 0);
        headerWriter.writeEnd(id, 0);

        size_t seqId = reader.getId(key);
        if (seqId == UINT_MAX) {
            Debug(Debug::ERROR) << "GFF entry not found in sequence database: " << name << "!\n";
            return EXIT_FAILURE;
        }

        ssize_t length = end - start;
        char *seq = reader.getData(seqId, 0);

        writer.writeStart(0);
        if (length > 0) {
            writer.writeAdd(seq + start, length + 1, 0);
        } else {
            for (size_t i = start; i >= start + length; i--) {
                revStr.append(1, Orf::complement(seq[i]));
            }
            writer.writeAdd(revStr.c_str(), revStr.size(), 0);
            revStr.clear();
        }
        writer.writeAdd("\n", 1, 0);
        writer.writeEnd(id, 0);

        entries_num++;
    }
    headerWriter.close();
    writer.close();
    headerReader.close();
    reader.close();
    file.close();

    return EXIT_SUCCESS;
}



