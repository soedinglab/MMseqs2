#include "Aggregation.h"
#include "Util.h"
#include "Debug.h"

#ifdef OPENMP
#include <omp.h>
#endif

Aggregation::Aggregation(const std::string &targetDbName, const std::string &resultDbName,
                         const std::string &outputDbName, unsigned int threads)
        : resultDbName(resultDbName), outputDbName(outputDbName), threads(threads) {
    std::string sizeDbName = targetDbName + "_member_to_set";
    std::string sizeDbIndex = targetDbName + "_member_to_set.index";
    targetSetReader = new DBReader<unsigned int>(sizeDbName.c_str(), sizeDbIndex.c_str());
    targetSetReader->open(DBReader<unsigned int>::NOSORT);
}

Aggregation::~Aggregation() {
    targetSetReader->close();
    delete targetSetReader;
}

// build a map with the value in [target column] field as a key and the rest of the line, cut in fields, as values
void Aggregation::buildMap(char *data, std::map<unsigned int, std::vector<std::vector<std::string>>> &dataToAggregate) {
    while (*data != '\0') {
        char *current = data;
        data = Util::skipLine(data);
        size_t length = data - current;
        std::string line(current, length - 1);
        if (line.empty() == true) {
            continue;
        }

        std::vector<std::string> columns = Util::split(line, "\t");
        unsigned int targetKey = (unsigned int) strtoull(columns[0].c_str(), NULL, 10);
        size_t setId = targetSetReader->getId(targetKey);
        if (setId == UINT_MAX) {
            Debug(Debug::ERROR) << "Invalid target database key " << columns[0] << ".\n";
            EXIT(EXIT_FAILURE);
        }
        char *data = targetSetReader->getData(setId);
        unsigned int setKey = (unsigned int) strtoull(data, NULL, 10);
        dataToAggregate[setKey].push_back(columns);
    }
}

int Aggregation::run() {
    std::string inputDBIndex = resultDbName + ".index";
    DBReader<unsigned int> reader(resultDbName.c_str(), inputDBIndex.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::string outputDBIndex = outputDbName + ".index";
    DBWriter writer(outputDbName.c_str(), outputDBIndex.c_str(), threads);
    writer.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(10 * 1024);

        std::map<unsigned int, std::vector<std::vector<std::string>>> dataToMerge;
#pragma omp for
        for (size_t i = 0; i < reader.getSize(); i++) {
            Debug::printProgress(i);
            dataToMerge.clear();

            unsigned int key = reader.getDbKey(i);
            buildMap(reader.getData(i), dataToMerge);
            prepareInput(key, thread_idx);
            
            for (std::map<unsigned int, std::vector<std::vector<std::string>>>::const_iterator it = dataToMerge.begin();
                 it != dataToMerge.end(); ++it) {
                unsigned int targetKey = it->first;
                std::vector<std::vector<std::string>> columns = it->second;
                buffer.append(aggregateEntry(columns, key, targetKey, thread_idx));
                buffer.append("\n");
            }
            writer.writeData(buffer.c_str(), buffer.length(), key, thread_idx);
            buffer.clear();
        }
    };
    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}
