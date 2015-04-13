#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <cstdlib>
#include <cstdio>

#include <string>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <algorithm>

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

struct blamap {
    int start;
    int end;
};

int findsorf(int argn, const char **argv)
{
    std::string usage;
    usage.append("Turns a cs219 FFindex database into the lagacy cs219 format.\n");
    usage.append("USAGE: <ffindexCS219InDB> <ffindexA3MInDB> <lagacyCS219Out>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");
    std::vector<MMseqsParameter> orf_par = {
        Parameters::PARAM_V
    };

    Parameters par;
    par.parseParameters(argn, argv, usage, orf_par, 3);

    Debug::setDebugLevel(par.verbosity);

    FILE* mapFile = fopen(par.db2.c_str(), "r");
    if(mapFile == NULL) { exit(EXIT_FAILURE); }

    /*std::string out_filename  = std::string(par.db3 + ".ffdata");
    std::string out_index_filename = std::string(par.db3 + ".ffindex");

    DBWriter writer(out_filename.c_str(), out_index_filename.c_str());
    writer.open();
*/
    std::map<int, int> speciesMap;
    std::string lastSpeciesName;
    int currentSpeciesIndex = 0;

    char * line = NULL;
    size_t length = 0;
    ssize_t read = 0;

    int start = 0, end = 0;

    while ((read = getline(&line, &length, mapFile)) != -1) {
        ffnchomp(line, read);
        std::vector<std::string> fields = Util::split(line, "\t");
        if(fields.size() < 2) {
            Debug(Debug::WARNING) << "Line does not contain enough entries!";
            continue;
        }

        speciesMap.emplace(atoi(fields[0].c_str()), currentSpeciesIndex);

        end = currentSpeciesIndex;
        if(fields[1] != lastSpeciesName) {
            lastSpeciesName = fields[1];
            start = currentSpeciesIndex++;
        }
    }
    fclose(mapFile);
    free(line);

    int totalSpeciesCount = currentSpeciesIndex;

    DBReader reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader::NOSORT);

    size_t num_clusters = reader.getSize();
    for (size_t i = 0; i < num_clusters; ++i) {
        const char* key = (const char*) reader.getDbKey(i);
        std::set<int> speciesSet;

        const char* clusters = (const char*)reader.getData(i);
        int entries = 0;
        std::istringstream iss(clusters);
        std::string line;
        while (std::getline(iss, line))
        {
            std::vector<std::string> clusterOrf = Util::split(line, "_");
            int speciesId = atoi(clusterOrf[0].c_str());
            speciesSet.insert(speciesMap[speciesId]);
            entries++;
        }

        if(entries <= 1)
            continue;

        float speciesInSet = (float)speciesSet.size();
        printf("%s\t%f \n", key, (speciesInSet/(float)totalSpeciesCount) * log2(speciesInSet / (float)totalSpeciesCount));
    }


    //writer.close();
    reader.close();

    return EXIT_SUCCESS;
}

