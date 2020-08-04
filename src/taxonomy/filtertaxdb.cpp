#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "TaxonomyExpression.h"

#ifdef OPENMP
#include <omp.h>
#endif

int filtertaxdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    NcbiTaxonomy* t = NcbiTaxonomy::openTaxonomy(par.db1);
    
    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    Debug::Progress progress(reader.getSize());
    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        char dbKey[255];
        TaxonomyExpression taxonomyExpression(par.taxonList, *t);

        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);

            writer.writeStart(thread_idx);
            while (*data != '\0') {
                Util::parseKey(data, dbKey);
                unsigned int taxon = Util::fast_atoi<unsigned int>(dbKey);
                if (taxonomyExpression.isAncestor(taxon)) {
                    char* nextData = Util::skipLine(data);
                    size_t dataSize = nextData - data;
                    writer.writeAdd(data, dataSize, thread_idx);
                }
                data = Util::skipLine(data);
            }
            writer.writeEnd(key, thread_idx);
        }
    }
    writer.close();
    reader.close();
    delete t;

    return EXIT_SUCCESS;
}
