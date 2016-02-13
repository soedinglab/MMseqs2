#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

#include <climits>
#include <algorithm>
#include <fstream>
#include <sstream>

//#ifdef OPENMP
//#include <omp.h>
//#endif

struct Edge {
    unsigned int identifier;
    unsigned int parent;

    Edge() {};
    Edge(unsigned int identifier, unsigned int parent) : identifier(identifier), parent(parent) {};

    static inline bool compareByIdentifier (const Edge& e1, const Edge& e2) {
        return e1.identifier < e2.identifier;
    }
};

std::vector<Edge>::const_iterator searchEdge(unsigned int identifier, const std::vector<Edge>& edges) {
    Edge search;
    search.identifier = identifier;
    std::vector<Edge>::const_iterator first = std::lower_bound(edges.begin(), edges.end(), search, Edge::compareByIdentifier);
 
    if (LIKELY(first != edges.end() && !(identifier < (*first).identifier))) {
        return first;
    } else {
        return edges.end();
    }
}

int clusteredges(int argc, const char *argv[]) {
    std::string usage("Extracts all child to parent edges from clustering N to clustering N+1.\n");
    usage.append("Written by Milot Mirdita <milot@mirdita.de>.\n\n");
    usage.append("USAGE: <edgeListOut> <clusteringSteps...> \n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.onlyverbosity, 1, true, true);
    Debug::setDebugLevel(par.verbosity);

    size_t clusterings = static_cast<size_t >(argc - 1);

    if (clusterings == 0) {
        Debug(Debug::ERROR) << "Need at least one cluster result to extract edges\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> *readers[clusterings];
    for (size_t i = 0; i < clusterings; ++i) {
        std::pair<std::string, std::string> name = Util::databaseNames(argv[i + 1]);
        readers[i] = new DBReader<unsigned int>(name.first.c_str(), name.second.c_str());
    }

//    if (FileUtil::fileExists(par.db1.c_str())) {
//        errno = EEXIST;
//        perror(par.db1.c_str());
//        EXIT(EXIT_FAILURE);
//    }
//
    std::vector<Edge> edges;

    for (size_t i = 0; i < clusterings; ++i) {
        DBReader<unsigned int> *reader = readers[i];
        reader->open(DBReader<unsigned int>::NOSORT);

        std::vector<Edge> tmp;

        for (size_t j = 0; j < reader->getSize(); ++j) {
            unsigned int parentKey = reader->getDbKey(j);

            // add the implicit root explicitly
            if (i == 0) {
                edges.emplace_back(parentKey, UINT_MAX);
            }

            std::istringstream ss(reader->getData(j));

            std::string line;
            std::getline(ss, line);
            while (std::getline(ss, line)) {
                if (line.length() < 1)
                    continue;
                
                unsigned int childKey = static_cast<unsigned int>(strtoul(line.c_str(), NULL, 10));
                edges.emplace_back(childKey, parentKey);
            }
        }

        reader->close();
    }

    for (size_t i = 0; i < clusterings; ++i) {
        delete readers[i];
    }

    std::sort(edges.begin(), edges.end(), Edge::compareByIdentifier);

    std::ofstream edgeList;
    edgeList.open(par.db1);

    if (edgeList.fail()) {
        Debug(Debug::ERROR) << "Could not open " << par.db1 << " for writing.";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db1 << "\n";
    for (std::vector<Edge>::const_iterator it = edges.begin(); it != edges.end(); ++it) {
        Edge edge = *it;
        if(edge.parent == UINT_MAX) {
            edgeList << edge.identifier << "\t{" << edge.identifier << "}\t0\n";
            continue;
        }

        std::vector<Edge> path;
        unsigned int depth = 0;
        unsigned int current = edge.parent;
        path.push_back(edge);

        while(current != UINT_MAX) {
            std::vector<Edge>::const_iterator edge = searchEdge(current, edges);
            if(edge != edges.end()) {
                current = (*edge).parent;
                path.push_back(*edge);
                depth++;
            } else {
                Debug(Debug::ERROR) << "Could not find parent " << current << "\n";
                EXIT(EXIT_FAILURE);
            }
        }

        std::ostringstream line;
        line << edge.identifier << "\t{";
        for(std::vector<Edge>::const_iterator p = path.begin(); p != path.end(); ++p) {
            Edge e = *p;
            line << e.identifier;
            if((p+1) != path.end())
                line << ",";
        }
        line << "}\t" << depth << "\n";

        edgeList << line.str();
    }

    edgeList.close();

    return EXIT_SUCCESS;
}
