//
// Created by lars on 02.06.15.
//

#include "convertfiles.h"
#include "Debug.h"
#include "DBWriter.h"
#include "Util.h"
#include <set>
#include <list>
#include <fstream>

convertfiles::convertfiles(std::string sequencedb, bool use_header) {
    if (use_header) {
        targetdb_header = new DBReader<unsigned int>(std::string(sequencedb + "_h").c_str(),
                                                     std::string(sequencedb + "_h.index").c_str());
        targetdb_header->open(DBReader<unsigned int>::NOSORT);
    } else {
        targetdb_header = NULL;
    }
    this->use_header = use_header;
}

convertfiles::~convertfiles() {
    if (targetdb_header != NULL) {
        targetdb_header->close();
        delete targetdb_header;
        targetdb_header = NULL;
    }
}

void convertfiles::convertFfindexToTsv(std::string clusteringfile, std::string suffix,
                                       std::string outputfile) {
    Debug(Debug::INFO) << "transfering: " << clusteringfile << " to " << outputfile << suffix << "cluster.tsv";

    std::string cluster_ffindex_indexfile = clusteringfile + ".index";

    DBReader<unsigned int> cluster_ffindex_reader(clusteringfile.c_str(),
                                                  cluster_ffindex_indexfile.c_str());
    cluster_ffindex_reader.open(DBReader<std::string>::SORT_BY_LENGTH);

    std::ofstream outfile_stream;
    outfile_stream.open(outputfile + suffix + "cluster.tsv");
    std::ofstream outfile_stream_clustersize;
    outfile_stream_clustersize.open(outputfile + suffix + "cluster_size.tsv");
    std::ofstream outfile_stream_cluster_summary;
    outfile_stream_cluster_summary.open(outputfile + suffix + "cluster_summary.tsv");

    outfile_stream_clustersize << "algorithm\tclusterid\tclustersize\n";
    outfile_stream << "algorithm\tclusterid\tid2\n";
    outfile_stream_cluster_summary << "algorithm\tclusternumber\tsingletons\n";
    int singletons = 0;
    char *idbuffer1 = new char[255 + 1];

    for (size_t i = 0; i < cluster_ffindex_reader.getSize(); ++i) {

        int clustersize = 0;
        char *data = cluster_ffindex_reader.getData(i);
        std::string idbuffer = "";
        std::string representative = "";
        while (*data != '\0') {
            Util::parseKey(data, idbuffer1);
            if (use_header) {
                idbuffer = getProteinNameForID(atoi(idbuffer1));
            } else {
                idbuffer = idbuffer1;
            }

            if (representative == "") {
                representative = idbuffer;
            }
            outfile_stream << suffix << "\t" << representative << "\t" << idbuffer << "\n";
            data = Util::skipLine(data);
            clustersize++;
        }
        if (clustersize == 1) {
            singletons++;
        }
        outfile_stream.flush();
        outfile_stream_clustersize << suffix << "\t" << representative << "\t" << clustersize << "\n";
    }
    outfile_stream_cluster_summary << suffix << "\t" << cluster_ffindex_reader.getSize() << "\t" << singletons << "\n";

    delete[] idbuffer1;

    cluster_ffindex_reader.close();

    outfile_stream.close();
    outfile_stream_clustersize.close();
    outfile_stream_cluster_summary.close();
}


void convertfiles::getAlignmentscoresForCluster(std::string clusteringfile, std::string alignmentfile,
                                                std::string outputfile) {
    //Debug(Debug::INFO) <<clusteringfile <<alignmentfile<<outputfile;

    std::string cluster_ffindex_indexfile = clusteringfile + ".index";
    std::string alignment_ffindex_indexfile = alignmentfile + ".index";

    DBReader<unsigned int> cluster_ffindex_reader(clusteringfile.c_str(),
                                                  cluster_ffindex_indexfile.c_str());
    cluster_ffindex_reader.open(DBReader<unsigned int>::SORT_BY_LENGTH);

    DBReader<unsigned int> alignment_ffindex_reader(alignmentfile.c_str(),
                                                    alignment_ffindex_indexfile.c_str());
    alignment_ffindex_reader.open(DBReader<unsigned int>::SORT_BY_LENGTH);

    std::ofstream outfile_stream;
    outfile_stream.open(outputfile);

    outfile_stream << "clusterid\tid2\taliscore\tqcov\tdbcov\tseqId\teval\n";

    const size_t LINE_BUFFER_SIZE = 255 + 1;
    char *idbuffer = new char[LINE_BUFFER_SIZE];
    char *linebuffer = new char[LINE_BUFFER_SIZE];

    for (size_t i = 0; i < cluster_ffindex_reader.getSize(); ++i) {


        unsigned int representative = cluster_ffindex_reader.getDbKey(i);
        char *data = cluster_ffindex_reader.getData(i);

        std::set<std::string> clusterset;


        while (*data != '\0') {
            Util::parseKey(data, idbuffer);
            clusterset.insert(std::string(idbuffer));
            data = Util::skipLine(data);
        }

        size_t dbKey = alignment_ffindex_reader.getId(representative);
        char *data_alignment = alignment_ffindex_reader.getData(dbKey);
        size_t data_alignment_length = alignment_ffindex_reader.getSeqLens(dbKey);
        if (data_alignment == NULL) {
            //  Debug(Debug::INFO) <<representative<<"\n";
            continue;
        }
        while (*data_alignment != '\0') {
            Util::parseKey(data_alignment, idbuffer);
            //Debug(Debug::INFO) <<idbuffer;
            if (clusterset.find(idbuffer) != clusterset.end()) {
                if (!Util::getLine(data_alignment, data_alignment_length, linebuffer, LINE_BUFFER_SIZE)) {
                    Debug(Debug::WARNING) << "Warning: Identifier was too long and was cut off!\n";
                }
                outfile_stream << getProteinNameForID(representative) << "\t" << linebuffer << "\n";
            }
            data_alignment = Util::skipLine(data_alignment);
        }
        outfile_stream.flush();

    }

    delete[] idbuffer;
    delete[] linebuffer;

    alignment_ffindex_reader.close();
    cluster_ffindex_reader.close();
    outfile_stream.close();

}

void convertfiles::convertDomainFileToFFindex(std::string domainscorefile, std::string domainIdentifierFile,
                                              std::string outputfile) {
    /*
     * Read in Domain identifiers
     * Line number (1 based index) is the index of the identifier in the domainfile
     */
    std::ifstream domainIdentifierFile_stream;
    domainIdentifierFile_stream.open(domainIdentifierFile);
    std::list<std::string> domainIdentifiers;
    std::string line;
    while (std::getline(domainIdentifierFile_stream, line)) {

        domainIdentifiers.push_back(std::string(line));

    }
    int number = 1;
    std::string *domainIdentifier = new std::string[domainIdentifiers.size() + 1];
    for (std::string id:domainIdentifiers) {
        domainIdentifier[number] = id;
        // Debug(Debug::INFO)<<id<<"\n";
        number++;
    }
    /*
     * read in scorefile , replace identifier and transfer to ffindex format
     */
    std::ifstream domainscorefile_stream;
    domainscorefile_stream.open(domainscorefile);

    std::string outputfileindex = outputfile + ".index";
    DBWriter dbw(outputfile.c_str(), outputfileindex.c_str());
    dbw.open();

    std::stringstream res;
    int lastid = -1;

    int a, b;
    double c;
    while (std::getline(domainscorefile_stream, line)) {
        std::istringstream iss(line);

        if (!(iss >> a >> b >> c)) {
            Debug(Debug::INFO) << line << "\n";
            break;
        }
        res << domainIdentifier[b] << "\t" << c << "\n";

        // domainIdentifier[a]<<"\t"<<domainIdentifier[b]<<"\t"<<res.str().length()<<"\n";

        if (lastid != -1 && a != lastid) {
            std::string cluResultsOutString = res.str();
            dbw.writeData(cluResultsOutString.c_str(), cluResultsOutString.length(), domainIdentifier[lastid].c_str());
            //clear res
            res.str("");
            res.clear();
        }
        lastid = a;
    }

    std::string cluResultsOutString = res.str();

    //save last entry
    dbw.writeData(cluResultsOutString.c_str(), cluResultsOutString.length(), domainIdentifier[a].c_str());
    dbw.close();

    delete[] domainIdentifier;
    domainIdentifierFile_stream.close();
}

void convertfiles::getDomainScoresForCluster(std::string clusteringfile, std::string alignmentfile,
                                             std::string outputfolder, std::string prefix, bool allagainstall,
                                             bool randomized) {
    // Debug(Debug::INFO) <<clusteringfile <<alignmentfile<< outputfolder;

    std::string cluster_ffindex_indexfile = clusteringfile + ".index";
    std::string alignment_ffindex_indexfile = alignmentfile + ".index";
    DBReader<unsigned int> cluster_ffindex_reader(clusteringfile.c_str(),
                                                   cluster_ffindex_indexfile.c_str());
    cluster_ffindex_reader.open(DBReader<std::string>::SORT_BY_LENGTH);
    DBReader<std::string> alignment_ffindex_reader(alignmentfile.c_str(),
                                                   alignment_ffindex_indexfile.c_str());
    alignment_ffindex_reader.open(DBReader<std::string>::SORT_BY_LENGTH);

    std::ofstream outfile_stream;
    outfile_stream.open(outputfolder + "/" + prefix + "domainscore.tsv");

    outfile_stream << "algorithm\tclusterid\tid2\tdomain_score\n";
    const size_t LINE_BUFFER_SIZE = 255 + 1;


    char *idbuffer1 = new char[LINE_BUFFER_SIZE];
    char *idbuffer2 = new char[LINE_BUFFER_SIZE];
    char *linebuffer = new char[LINE_BUFFER_SIZE];

    for (size_t i = 0; i < cluster_ffindex_reader.getSize(); ++i) {
        std::string representative = getProteinNameForID(cluster_ffindex_reader.getDbKey(i));
        char *data = cluster_ffindex_reader.getData(i);

        std::set<std::string> clusterset;
        std::set<std::string> clusterset2;

        while (*data != '\0') {
            Util::parseKey(data, idbuffer1);
            std::string idbuffer = getProteinNameForID(atoi(idbuffer1));
            clusterset.insert(std::string(idbuffer));
            data = Util::skipLine(data);
        }

        size_t dbKey = alignment_ffindex_reader.getId(representative);
        char *data_alignment = alignment_ffindex_reader.getData(dbKey);
        size_t data_alignment_length = alignment_ffindex_reader.getSeqLens(dbKey);
        if (data_alignment == NULL) {
            //  Debug(Debug::INFO) <<representative<<"\n";
            continue;
        }
        // Debug(Debug::INFO) <<representative<<"\n";
        while (*data_alignment != '\0') {
            Util::parseKey(data_alignment, idbuffer2);
            // Debug(Debug::INFO) <<idbuffer;
            if (clusterset.find(idbuffer2) != clusterset.end()) {
                clusterset2.insert(std::string(idbuffer2));
                if (!Util::getLine(data_alignment, data_alignment_length, linebuffer, LINE_BUFFER_SIZE)) {
                    Debug(Debug::WARNING) << "Warning: Identifier was too long and was cut off!\n";
                }
                outfile_stream << prefix << "\t" << representative << "\t" << linebuffer << "\n";
            }
            data_alignment = Util::skipLine(data_alignment);
        }
        outfile_stream.flush();


        for (std::string id :clusterset) {
            strncpy(idbuffer2, id.c_str(), id.length());
            idbuffer2[id.length()] = '\0';
            if (strcmp(representative.c_str(), id.c_str()) == 0 || clusterset2.find(id) != clusterset2.end()) {

            } else {
                if (alignment_ffindex_reader.getDataByDBKey(idbuffer2) != NULL) {
                    outfile_stream << prefix << "\t" << representative << "\t" << id << "\t" << "0" << "\n";
                }
            }

        }
        outfile_stream.flush();
    }

    delete[] idbuffer1;
    delete[] idbuffer2;
    delete[] linebuffer;

    alignment_ffindex_reader.close();
    cluster_ffindex_reader.close();
    outfile_stream.close();
}


std::string convertfiles::getProteinNameForID(unsigned int dbKey) {
    char *header_data = targetdb_header->getDataByDBKey(dbKey);
    return Util::parseFastaHeader(header_data);
}
