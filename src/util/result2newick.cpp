#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

#include <climits>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <deque>
#include <set>

struct Edge {
    unsigned int identifier;
    unsigned int parent;
    unsigned int step;

    Edge() {};
    Edge(unsigned int identifier, unsigned int parent, unsigned int step) : identifier(identifier), parent(parent), step(step) {};

    static inline bool compareByIdentifier (const Edge& e1, const Edge& e2) {
        return e1.identifier < e2.identifier;
    }
};

class Node {
public:
    unsigned int identifier;
    unsigned int step;

    std::set<Node> children;

    explicit Node(unsigned int identifier, unsigned int step) : identifier(identifier), step(step) {}

    bool operator<( const Node& other ) const
    {
        return identifier < other.identifier;
    }
};

std::string nodeAsNewick(const Node* node, int depth)
{
    std::ostringstream ss;
    if (node->children.size() == 0) {
        ss << node->identifier << ":" << node->step;
    } else {
        ss << "(";
        for (std::set<Node>::const_iterator it = node->children.begin(); it != node->children.end(); ++it) {
            const Node* n = &(*it);
            ss << nodeAsNewick(n, depth+1);
            if(std::distance(it, node->children.end()) != 1) {
               ss << ",";
            }
        }

        ss << ")";

        if (depth > 0)
            ss << node->identifier << ":" << node->step;
        else
            ss << node->identifier;
    }

    return ss.str();
}

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

int result2newick(int argc, const char *argv[]) {
    std::string usage("Extracts clustering relationship from clustering steps into Newick trees.\n");
    usage.append("Written by Milot Mirdita <milot@mirdita.de>.\n\n");
    usage.append("USAGE: <edgeListOutDB> <clusteringSteps0...n> \n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.onlyverbosity, 1, true, true);

    size_t clusterings = static_cast<size_t>(argc - 1);

    if (clusterings == 0) {
        Debug(Debug::ERROR) << "Need at least one cluster result to extract edges\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> **readers = new DBReader<unsigned int> *[clusterings];
    for (size_t i = 0; i < clusterings; ++i) {
        std::pair<std::string, std::string> name = Util::databaseNames(argv[i + 1]);
        readers[i] = new DBReader<unsigned int>(name.first.c_str(), name.second.c_str());
    }

    std::vector<Edge> edges;
    for (size_t i = 0; i < clusterings; ++i) {
        DBReader<unsigned int>& reader = *readers[i];
        reader.open(DBReader<unsigned int>::NOSORT);

        std::vector<Edge> tmp;

        for (size_t j = 0; j < reader.getSize(); ++j) {
            unsigned int parentKey = reader.getDbKey(j);

            std::istringstream ss(reader.getData(j));

            std::string line;
            std::getline(ss, line);
            while (std::getline(ss, line)) {
                if (line.length() < 1)
                    continue;
                
                unsigned int childKey = static_cast<unsigned int>(strtoul(line.c_str(), NULL, 10));
                edges.emplace_back(childKey, parentKey, i);
            }
        }

        reader.close();
    }

    std::sort(edges.begin(), edges.end(), Edge::compareByIdentifier);

    
    std::vector<std::deque<std::pair<unsigned int, unsigned int>>> paths;

    DBReader<unsigned int>& rootReader = *readers[0];
    rootReader.open(DBReader<unsigned int>::NOSORT);

    for (size_t i = 0; i < rootReader.getSize(); ++i) {
        unsigned int id = rootReader.getDbKey(i);

        std::deque<std::pair<unsigned int, unsigned int>> path;
        path.emplace_front(id, 0);

        unsigned int step = 0;
        std::vector<Edge>::const_iterator edge;
        while((edge = searchEdge(id, edges)) != edges.end()) {
            Edge e = *edge;
            id = e.parent;
//            step = e.step;
            unsigned  int skip = e.step - step;
            while(skip > 0) {
                path.emplace_front(id, e.step);
                skip--;
            }


            step = e.step;
        }

        for(unsigned int j = step; j < clusterings - 1; ++j) {
            path.emplace_front(id, j + 1);
        }

        paths.push_back(path);
    }
    edges.clear();

    rootReader.close();

    for (size_t i = 0; i < clusterings; ++i) {
        delete readers[i];
    }

    std::set<Node> root;
    for(std::vector<std::deque<std::pair<unsigned int, unsigned int>>>::const_iterator i = paths.begin(); i != paths.end(); ++i) {
        std::set<Node>* current = &root;
        for(std::deque<std::pair<unsigned int, unsigned int>>::const_iterator j = (*i).begin(); j != (*i).end(); ++j) {
            Node search((*j).first, (*j).second);
            std::pair<std::set<Node>::iterator, bool> result = current->insert(search);
            current = const_cast<std::set<Node>*>(&((result.first)->children));
        }
    }

    paths.clear();

    DBWriter writer(par.db1.c_str(), par.db1Index.c_str(), 1);
    writer.open();

    for (std::set<Node>::const_iterator it = root.begin(); it != root.end(); ++it) {
        Node node = *it;
        std::string newick = nodeAsNewick(&node, 0);

        std::set<unsigned int> identifiers;
        std::deque<Node> visit;
        visit.push_back(node);
        while(visit.size() > 0) {
            Node n = visit.front();
            visit.pop_front();
            identifiers.insert(n.identifier);
            visit.insert(visit.end(), n.children.begin(), n.children.end());
        }

        for (std::set<unsigned int>::const_iterator it = identifiers.begin(); it != identifiers.end(); ++it) {
            writer.writeData(newick.c_str(), newick.length(), SSTR(*it).c_str());
        }
    }

    writer.close();
    delete [] readers;
    return EXIT_SUCCESS;
}
