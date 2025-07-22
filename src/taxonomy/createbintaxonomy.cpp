#include "Debug.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "NcbiTaxonomy.h"

int createbintaxonomy(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, 0, 0);
    NcbiTaxonomy taxonomy(par.db1, par.db2, par.db3);
    std::pair<char*, size_t> serialized = NcbiTaxonomy::serialize(taxonomy);
    FILE* handle = fopen(par.db4.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << par.db4 << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    fwrite(serialized.first, serialized.second, sizeof(char), handle);
    fclose(handle);
    free(serialized.first);
    EXIT(EXIT_SUCCESS);
}

int createdmptaxonomy(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, 0, 0);
    NcbiTaxonomy* taxonomy = NcbiTaxonomy::openTaxonomy(par.db1);


    std::string nodesPath = par.db2 + "_nodes.dmp";
    FILE *nodesFp = FileUtil::openAndDelete(nodesPath.c_str(), "w");
    if (!nodesFp) {
        Debug(Debug::ERROR) << "Could not open " << nodesPath << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    std::string namesPath = par.db2 + "_names.dmp";
    FILE *namesFp = FileUtil::openAndDelete(namesPath.c_str(), "w");
    if (!namesFp) {
        Debug(Debug::ERROR) << "Could not open " << namesPath << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    for (size_t i = 0; i < taxonomy->maxNodes; ++i) {
        const TaxonNode &tn = taxonomy->taxonNodes[i];
        int res = fprintf(
            nodesFp,
            "%d\t|\t%d\t|\t%s\t|\t\n",
            tn.taxId,
            tn.parentTaxId,
            taxonomy->getString(tn.rankIdx)
        );
        if (res < 0) {
            Debug(Debug::ERROR) << "Write error while writing " << nodesPath << "\n";
            EXIT(EXIT_FAILURE);
        }

        const char *name = taxonomy->getString(tn.nameIdx);
        res = fprintf(
            namesFp,
            "%d\t|\t%s\t|\t\t|\tscientific name\t|\n",
            tn.taxId,
            name
        );
        if (res < 0) {
            Debug(Debug::ERROR) << "Write error while writing " << nodesPath << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    if (fclose(nodesFp) != 0) {
        Debug(Debug::ERROR) << "Could not close " << nodesPath << "\n";
        fclose(namesFp);
        EXIT(EXIT_FAILURE);
    }
    if (fclose(namesFp) != 0) {
        Debug(Debug::ERROR) << "Could not close " << namesPath << "\n";
        EXIT(EXIT_FAILURE);
    }

    EXIT(EXIT_SUCCESS);
}
