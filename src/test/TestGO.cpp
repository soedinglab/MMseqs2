//
// Created by lars on 12.04.15.
//

#include "Util.h"
#include "DistanceCalculator.h"
#include "convertfiles.h"
#include "DBWriter.h"
#include "Parameters.h"
#include "CompareGOTerms.h"
#include "Debug.h"

#include <cmath>
#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif

void printHelp();
std::string getProteinNameForID(DBReader<unsigned int>* targetdb_header,unsigned int dbKey);

int main(int argc, char **argv)
{

    if (argc < 1) {
        printHelp();

    }
    for(int i=0;i<argc;i++){
        Debug(Debug::INFO)<<argv[i]<<" ";
    }
    Debug(Debug::INFO)<<"\n";
   // for (int k = 0; k < 8; ++k) {
    //    Debug(Debug::INFO) << argv[k] <<"\n";
    //}
    //////////
    /////GO-Evaluation
    /////////
    if(strcmp(argv[1],"-go")==0){
        Debug(Debug::INFO) <<"GO-Evaluation" <<"\n";

        if(argc != 10){
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }

        std::string gofolder=argv[2];
        std::string uniprot_go_folder=argv[3];
        std::string clustering_file=argv[4];
        std::string prefix=argv[5];
        std::string outputfolder=argv[6];
        char * comparisonmode=argv[7];
        char *randommode =argv[8];
        std::string sequencedb=argv[9];

        bool allagainstall=false;
        bool randomized=false;

        if(strcmp(comparisonmode,"yes")==0){
            allagainstall=true;
            Debug(Debug::INFO)<<"all against all comparison";
        }
        if(strcmp(randommode,"yes")==0){
            randomized=true;
            Debug(Debug::INFO)<<"randomized representative comparison";
        }
        bool use_header=true;
        if(sequencedb=="no"){
            use_header=false;
            Debug(Debug::INFO)<<"do not use sequencedb for mapping, uniref mode";
        }
        //"-go <gofolder> <prot_go_folder> <clustering_file> <prefix> <outputfolder>"
        //std::string gofolder="/home/lars/masterarbeit/data/GO/db/";
        //std::string uniprot_go_folder="/home/lars/masterarbeit/data/uniprot/release-2015_04/uniprot_go/";
        //"/home/lars/masterarbeit/data/uniprot/release-2015_04/evaluation"
          //      "/home/lars/masterarbeit/db/sprot/uniprot_sprot_s4_affinity"

        std::string* goCategories= new std::string[3];
        goCategories[0]="_C";
        goCategories[1]="_F";
        goCategories[2]="_P";
        std::string* evidenceCategories= new std::string[3];
        evidenceCategories[0]="";
        evidenceCategories[1]="_EXP";
        evidenceCategories[2]="_NON-IEA";

        for (int j = 0; j <3 ; ++j) {
            for (int i = 0; i < 3; ++i) {
                CompareGOTerms *go = new CompareGOTerms(gofolder + "go-fasta_db" + goCategories[i],
                                                        gofolder + "go-fasta_db" + goCategories[i] + ".index",
                                                        uniprot_go_folder + "uniprot_sprot.dat_go_db" +
                                                        evidenceCategories[j] + goCategories[i],
                                                        uniprot_go_folder + "uniprot_sprot.dat_go_db" +
                                                        evidenceCategories[j] + goCategories[i] + ".index",
                                                        outputfolder,sequencedb,use_header);
                go->init();
                //go->all_against_all_comparison();
                //  go->all_against_all_comparison_proteinset();
                go->run_evaluation_mmseqsclustering(clustering_file,
                                                    clustering_file+".index",
                                                    prefix, evidenceCategories[j] + goCategories[i],allagainstall,randomized);
                go->~CompareGOTerms();
            }
        }
    }
    //////////
    /////Protein Name
    /////////
    else if (strcmp(argv[1],"-pn")==0){
        if(argc != 9){
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string protname_db=argv[2];
        std::string protname_db_indexfile=protname_db+".index";
        std::string cluster_ffindex=argv[3];
        std::string cluster_ffindex_indexfile=cluster_ffindex+".index";
        std::string fileprefix=argv[4];
        std::string evaluationfolder=argv[5];
        char * comparisonmode=argv[6];
        char * randomnmode=argv[7];
        std::string sequencedb=argv[8];
        bool allagainstall=false;
        bool randomized=false;
        if(strcmp(comparisonmode,"yes")==0){
            allagainstall=true;
            Debug(Debug::INFO)<<"all against all comparison";
        }
        if(strcmp(randomnmode,"yes")==0){
            randomized=true;
            Debug(Debug::INFO)<<"randomized representative comparison";
        }
        bool use_header=true;
        if(sequencedb=="no"){
            use_header=false;
            Debug(Debug::INFO)<<"do not use sequencedb for mapping, uniref mode";
        }
            Debug(Debug::INFO) << "running protein name evaluation";

        Debug(Debug::INFO) << "Opening clustering database...\n";
        DBReader<unsigned int>* cluster_ffindex_reader = new DBReader<unsigned int>(cluster_ffindex.c_str(), cluster_ffindex_indexfile.c_str());
        cluster_ffindex_reader->open(DBReader<std::string>::SORT_BY_LENGTH);
        Debug(Debug::INFO) << "Opening clustering database...\n";
        DBReader<std::string>* protname_db_reader = new DBReader<std::string>(protname_db.c_str(), protname_db_indexfile.c_str());
        protname_db_reader->open(DBReader<std::string>::NOSORT);
        DBReader<unsigned int> *targetdb_header=NULL;
        if(use_header) {
            targetdb_header = new DBReader<unsigned int>(std::string(sequencedb + "_h").c_str(),
                                                                                 std::string(sequencedb +
                                                                                             "_h.index").c_str());
            targetdb_header->open(DBReader<std::string>::NOSORT);
        }
        //files
        std::ofstream clusters_full_file;
        clusters_full_file.open(evaluationfolder+fileprefix+"clusters_allproteinnamescores.go");
        clusters_full_file << "algorithm\tclusterid\tid1\tid2\tproteinnamescore\n";



        for (size_t i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


                std::string representative="";

            char *data = cluster_ffindex_reader->getData(i);
            char *idbuffer1 = new char[255 + 1];
            double sumofscore =0;
            double minscore=1;
            double maxscore=0;

            std::list<std::string> idswithproteinname;


            //Debug(Debug::INFO) << representative << "\t" << "not available" << "\n";

            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                std::string idbuffer ="";
                if(use_header) {
                   idbuffer = getProteinNameForID(targetdb_header, atoi(idbuffer1));
                }else{
                     idbuffer =idbuffer1;
                }

                if(representative==""){
                    representative=idbuffer;
                }
                if (protname_db_reader->getDataByDBKey(idbuffer) != NULL) {
                    idswithproteinname.push_back(std::string(idbuffer));

                } else {
                    //Debug(Debug::INFO) << representative << "\t" << idbuffer << "\t" << "not available" <<"\n";
                }

                data = Util::skipLine(data);
            }

            for(std::string id1:idswithproteinname){
                //Debug(Debug::INFO) <<id1<<"\t";
                if(randomized){
                    std::list<std::string>::iterator i = idswithproteinname.begin();
                    std::advance(i, rand()% idswithproteinname.size());
                    id1=*i;
                }

                //Debug(Debug::INFO) <<id1<<"\n";
                for(std::string id2:idswithproteinname){
                    if (std::string(id1) != std::string(id2)) {
                        char*proteinName1 =protname_db_reader->getDataByDBKey((char *) id1.c_str());
                        char *proteinName2 =protname_db_reader->getDataByDBKey((char *) id2.c_str());
                        double levenshteinScore = DistanceCalculator::levenshteinDistance(proteinName1, proteinName2);
                        levenshteinScore = 1- (levenshteinScore /std::max(strlen(proteinName1),strlen(proteinName2)));
                        sumofscore += levenshteinScore;
                        minscore=std::min(levenshteinScore,minscore);
                        maxscore=std::max(levenshteinScore,maxscore);
                        clusters_full_file << fileprefix << "\t"<<  representative.c_str() << "\t" << id1.c_str() << "\t" << id2.c_str() << "\t" <<
                                                                                                                    levenshteinScore << "\n";
                    }
                }
                if(!allagainstall){
                    break;
                }

            }

        }




        clusters_full_file.close();

        cluster_ffindex_reader->close();
        cluster_ffindex_reader->~DBReader();
        protname_db_reader->~DBReader();



    }
        //////////
        /////Key-WordEvaluation
        /////////
    else if (strcmp(argv[1],"-kw")==0){
        if(argc != 9){
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string keyword_db=argv[2];
        std::string keyword_indexfile=keyword_db+".index";
        std::string cluster_ffindex=argv[3];
        std::string cluster_ffindex_indexfile=cluster_ffindex+".index";
        std::string fileprefix=argv[4];
        std::string evaluationfolder=argv[5];
        char * comparisonmode=argv[6];
        char * randomnmode=argv[7];
        std::string sequencedb=argv[8];
        bool allagainstall=false;
        bool randomized=false;
        if(strcmp(comparisonmode,"yes")==0){
            allagainstall=true;
            Debug(Debug::INFO)<<"all against all comparison";
        }
        if(strcmp(randomnmode,"yes")==0){
            randomized=true;
            Debug(Debug::INFO)<<"randomized representative comparison";
        }
        bool use_header=true;
        if(sequencedb=="no"){
            use_header=false;
            Debug(Debug::INFO)<<"do not use sequencedb for mapping, uniref mode";
        }
        Debug(Debug::INFO) << "running keyword evaluation";

        Debug(Debug::INFO) << "Opening clustering database...\n";
        DBReader<unsigned int>* cluster_ffindex_reader = new DBReader<unsigned int>(cluster_ffindex.c_str(), cluster_ffindex_indexfile.c_str());
        cluster_ffindex_reader->open(DBReader<std::string>::NOSORT);
        Debug(Debug::INFO) << "Opening clustering database...\n";
        DBReader<std::string>* protname_db_reader = new DBReader<std::string>(keyword_db.c_str(), keyword_indexfile.c_str());
        protname_db_reader->open(DBReader<std::string>::NOSORT);
        DBReader<unsigned int> *targetdb_header=NULL;
        if(use_header) {
            targetdb_header = new DBReader<unsigned int>(std::string(sequencedb + "_h").c_str(),
                                                                                 std::string(sequencedb +
                                                                                             "_h.index").c_str());
            targetdb_header->open(DBReader<std::string>::NOSORT);
        }

        //files
        std::ofstream clusters_full_file;
        clusters_full_file.open(evaluationfolder+fileprefix+"clusters_keywordscores.go");
        clusters_full_file << "algorithm\tclusterid\tid1\tid2\tkeywordscore\n";



        for (unsigned int i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


            std::string representative="";
            char *data = cluster_ffindex_reader->getData(i);
            char *idbuffer1 = new char[255 + 1];
            double sumofscore =0;
            double minscore=1;
            double maxscore=0;

            std::list<std::string> idswithkeyword;


            //Debug(Debug::INFO) << representative << "\t" << "not available" << "\n";

            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                std::string idbuffer ="";
                if(use_header) {
                    idbuffer = getProteinNameForID(targetdb_header, atoi(idbuffer1));
                }else{
                    idbuffer =idbuffer1;
                }

                if(representative==""){
                    representative=idbuffer;
                }

                if (protname_db_reader->getDataByDBKey(idbuffer) != NULL) {
                    idswithkeyword.push_back(std::string(idbuffer));


                } else {
                    //Debug(Debug::INFO) << representative << "\t" << idbuffer << "\t" << "not available" <<"\n";
                }

                data = Util::skipLine(data);
            }
            for(std::string id1:idswithkeyword){
                if(randomized){
                    std::list<std::string>::iterator i = idswithkeyword.begin();
                    std::advance(i, rand()% idswithkeyword.size());
                    id1=*i;
                }
                for(std::string id2:idswithkeyword){
                    if (std::string(id1) != std::string(id2)) {
                        char* seq1=protname_db_reader->getDataByDBKey((char *) id1.c_str());
                        char * seq2=protname_db_reader->getDataByDBKey((char *) id2.c_str());
                        double score = DistanceCalculator::keywordDistance( seq1,  seq2);
                        sumofscore += score;
                        minscore=std::min(score,minscore);
                        maxscore=std::max(score,maxscore);
                        clusters_full_file << fileprefix << "\t"<<  representative << "\t" << id1 << "\t" << id2 << "\t" << score << "\n";
                    }
                }
                if(!allagainstall){
                    break;
                }

            }
            /*double averagescore;
            if(allagainstall){
                averagescore=(sumofscore / (idswithgo.size()*idswithgo.size()-idswithgo.size()));
            }else{
                averagescore=sumofscore /(idswithgo.size()-1);
            }*/


        }




        clusters_full_file.close();

        cluster_ffindex_reader->close();
        cluster_ffindex_reader->~DBReader();
        protname_db_reader->~DBReader();



    }else if (strcmp(argv[1],"-cs")==0) {
        if (argc != 6) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string clusteringfile = argv[2];
        std::string alignmentfile = argv[3];
        std::string outputfile = argv[4];
        std::string sequencedb=argv[5];

        convertfiles *cf = new convertfiles(sequencedb,true);
        cf->getAlignmentscoresForCluster(clusteringfile,alignmentfile,outputfile);

    }else if (strcmp(argv[1],"-df")==0) {
        if (argc != 6) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string domainscorefile = argv[2];
        std::string domainIdentifierFile = argv[3];
        std::string outputfile = argv[4];
        std::string sequencedb=argv[5];

        convertfiles *cf = new convertfiles(sequencedb,true);
        cf->convertDomainFileToFFindex(domainscorefile,domainIdentifierFile,outputfile);

    }else if (strcmp(argv[1],"-ds")==0) {
        if (argc != 9) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string clusteringfile = argv[2];
        std::string alignmentfile = argv[3];
        std::string suffix = argv[4];
        std::string outputfile = argv[5];
        char * comparisonmode=argv[6];
        char * randomnmode=argv[7];
        std::string sequencedb=argv[8];
        bool allagainstall=false;
        bool randomized=false;
        if(strcmp(comparisonmode,"yes")==0){
            allagainstall=true;
            Debug(Debug::INFO)<<"all against all comparison";
        }
        if(strcmp(randomnmode,"yes")==0){
            randomized=true;
            Debug(Debug::INFO)<<"randomized representative comparison";
        }
        convertfiles *cf = new convertfiles(sequencedb,true);
        cf->getDomainScoresForCluster(clusteringfile,alignmentfile,outputfile,suffix,allagainstall,randomized);

    }else if (strcmp(argv[1],"-clusterToTsv")==0) {
        if (argc != 6) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string clusteringfile = argv[2];
        std::string suffix = argv[3];
        std::string outputfolder = argv[4];
        std::string sequencedb=argv[5];

        bool use_header=true;
        if(sequencedb=="no"){
            use_header=false;
            Debug(Debug::INFO)<<"do not use sequencedb for mapping, uniref mode";
        }
        convertfiles *cf = new convertfiles(sequencedb,use_header);
        cf->convertFfindexToTsv(clusteringfile, suffix, outputfolder);

    }else if(strcmp(argv[1],"-af")==0){
        if (argc != 7) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        for(int i=0;i<argc;i++){
            Debug(Debug::INFO)<<argv[i]<<" ";
        }
        Debug(Debug::INFO)<<"\n";
        std::string seqDbfile = argv[2];
        std::string alignmentfile = argv[3];
        std::string outputfile = argv[4];
        float threshold = atof(argv[5]);
        int similarityScoreType = atof(argv[6]);
int numberofthreads=1;
#ifdef OPENMP
        omp_set_num_threads(numberofthreads);
#endif
        DBReader<unsigned int>* seqDbr=new DBReader<unsigned int>(seqDbfile.c_str(),(seqDbfile+".index").c_str());
        seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);
        DBReader<unsigned int>* alnDbr=new DBReader<unsigned int>(alignmentfile.c_str(),(alignmentfile+".index").c_str());
        alnDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);
        DBWriter* dbw = new DBWriter(outputfile.c_str(), (outputfile+".index").c_str(),numberofthreads);
        dbw->open();

        const size_t BUFFER_SIZE = 1000000;
        const size_t LINE_BUFFER_SIZE = 255+1;

#pragma omp parallel
        {
            char* outBuffer = new char[BUFFER_SIZE];
            char *idbuffer1 = new char[LINE_BUFFER_SIZE];
            char *linebuffer = new char[LINE_BUFFER_SIZE];
            char *similarity = new char[LINE_BUFFER_SIZE];

#pragma omp for schedule(dynamic, 100)
            for (unsigned int i = 0; i < seqDbr->getSize(); i++) {
                Debug::printProgress(i);
                int thread_idx = 0;
#ifdef OPENMP
                thread_idx = omp_get_thread_num();
#endif
                // for (int i = 0; i < seqDbr->getSize(); i++) {
                char *data = alnDbr->getData(i);
                size_t dataLength = alnDbr->getSeqLens(i);
                std::string cluResultsOutString = std::string("");
                while (*data != '\0') {
                    Util::parseKey(data, idbuffer1);
                    if(!Util::getLine(data, dataLength, linebuffer, LINE_BUFFER_SIZE)) {
                        Debug(Debug::WARNING) << "Warning: Identifier was too long and was cut off!\n";
                    }
                    //get similarityscore
                    float factor = 1;
                    float similarityscore=0.0;
                    if (similarityScoreType == Parameters::APC_SEQID) {
                        Util::parseByColumnNumber(data, similarity, 2); //column 4 = sequence identity
                        similarityscore = atof(std::string(similarity).c_str()) * factor;
                    }

                    if (similarityscore < threshold) {
                        data = Util::skipLine(data);
                        continue;
                    }
                    cluResultsOutString = cluResultsOutString + linebuffer + "\n";
                    data = Util::skipLine(data);
                }

                // Debug(Debug::INFO)<<data;

                const char *cluResultsOutData = cluResultsOutString.c_str();
                if (BUFFER_SIZE < strlen(cluResultsOutData)) {
                    Debug(Debug::ERROR) << "Tried to process the alignment list for the query " << alnDbr->getDbKey(i)
                    << " , length of the list = " << "\n";
                    Debug(Debug::ERROR) << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " <<
                    cluResultsOutString.length()
                    << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
                    continue;
                }
                memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length() * sizeof(char));
                dbw->writeData(outBuffer, cluResultsOutString.length(), std::to_string(alnDbr->getDbKey(i)).c_str(),
                               thread_idx);

                //  data = Util::skipLine(data);
                //  }
            }
        }
        alnDbr->close();
        seqDbr->close();
        dbw->close();
        alnDbr->~DBReader();
        seqDbr->~DBReader();
        dbw->~DBWriter();
    }else{
        printHelp();
        Debug(Debug::INFO)<<DistanceCalculator::levenshteinDistance("bla","bla21");
    }




    }

void printHelp() {
    std::string usage("\nEvaluation commands\n");
    usage.append("-go <gofolder> <prot_go_folder> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) ><yes : randomized representative choice |no : representative against all(default) > \n");
    usage.append("-pn <prot_name_db> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) ><yes : randomized representative choice |no : representative against all(default) > \n");
    usage.append("-kw <keyword_db> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) ><yes : randomized representative choice |no : representative against all(default) > \n");
    usage.append("-cs <clustering_file> <alignment_file> <outputfile> \n");
    usage.append("-df <domainscorefile> <domainIdentifierFile> <outputfile> \n");
    usage.append("-ds <clustering_file> <domainscorefile> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) ><yes : randomized representative choice |no : representative against all(default) > \n");
    usage.append("-clusterToTsv <clustering_file> <prefix> <outputfolder> <sequencedb>\n");
    usage.append("-af  <seqDbfile> <alignmentfile> <outputfile> <cutoff> <scoretype: 1-5>\n");

    Debug(Debug::INFO) << usage << "\n";
    EXIT(EXIT_FAILURE);
}


std::string getProteinNameForID(DBReader<unsigned int>* targetdb_header, unsigned int dbKey){
    char * header_data = targetdb_header->getDataByDBKey(dbKey);
    std::string parsedDbkey = Util::parseFastaHeader(header_data);
    return parsedDbkey;
}
