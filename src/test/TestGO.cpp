//
// Created by lars on 12.04.15.
//

#include <Util.h>
#include <DistanceCalculator.h>
#include <convertfiles.h>
#include <DBWriter.h>
#include <Parameters.h>
#include "CompareGOTerms.h"


void printHelp();

int main(int argc, char **argv)
{

    if (argc < 1) {
        printHelp();

    }
   // for (int k = 0; k < 8; ++k) {
    //    Debug(Debug::INFO) << argv[k] <<"\n";
    //}
    //////////
    /////GO-Evaluation
    /////////
    if(strcmp(argv[1],"-go")==0){
        Debug(Debug::INFO) <<"GO-Evaluation" <<"\n";

        if(argc != 8){
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }

        std::string gofolder=argv[2];
        std::string uniprot_go_folder=argv[3];
        std::string clustering_file=argv[4];
        std::string prefix=argv[5];
        std::string outputfolder=argv[6];
        char * comparisonmode=argv[7];
        bool allagainstall=false;
        if(strcmp(comparisonmode,"yes")==0){
            allagainstall=true;
            Debug(Debug::INFO)<<"all against all comparison";
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
                                                        outputfolder);
                go->init();
                //go->all_against_all_comparison();
                //  go->all_against_all_comparison_proteinset();
                go->run_evaluation_mmseqsclustering(clustering_file,
                                                    clustering_file+".index",
                                                    prefix, evidenceCategories[j] + goCategories[i],allagainstall);
                go->~CompareGOTerms();
            }
        }
    }
    //////////
    /////Protein Name
    /////////
    else if (strcmp(argv[1],"-pn")==0){
        if(argc != 7){
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
        bool allagainstall=false;
        if(strcmp(comparisonmode,"yes")==0){
            allagainstall=true;
            Debug(Debug::INFO)<<"all against all comparison";
        }
            Debug(Debug::INFO) << "running protein name evaluation";

        Debug(Debug::INFO) << "Opening clustering database...\n";
        DBReader* cluster_ffindex_reader = new DBReader(cluster_ffindex.c_str(), cluster_ffindex_indexfile.c_str());
        cluster_ffindex_reader->open(DBReader::SORT);
        Debug(Debug::INFO) << "Opening clustering database...\n";
        DBReader* protname_db_reader = new DBReader(protname_db.c_str(), protname_db_indexfile.c_str());
        protname_db_reader->open(DBReader::NOSORT);


        //files
        std::ofstream clusters_full_file;
        clusters_full_file.open(evaluationfolder+fileprefix+"clusters_allproteinnamescores.go");
        clusters_full_file << "algorithm\tclusterid\tid1\tid2\tproteinnamescore\n";



        for (int i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


            char *representative=cluster_ffindex_reader->getDbKey(i);
            char *data = cluster_ffindex_reader->getData(i);
            char *idbuffer = new char[255 + 1];
            double sumofscore =0;
            double minscore=1;
            double maxscore=0;

            std::list<std::string> idswithgo;


            //Debug(Debug::INFO) << representative << "\t" << "not available" << "\n";

            while (*data != '\0') {
                Util::parseKey(data, idbuffer);

                if (protname_db_reader->getDataByDBKey(idbuffer) != NULL) {
                    idswithgo.push_back(std::string(idbuffer));

                } else {
                    //Debug(Debug::INFO) << representative << "\t" << idbuffer << "\t" << "not available" <<"\n";
                }

                data = Util::skipLine(data);
            }
            for(std::string id1:idswithgo){
                for(std::string id2:idswithgo){
                    if (std::string(id1) != std::string(id2)) {
                        char* seq1=protname_db_reader->getDataByDBKey((char *) id1.c_str());
                        char * seq2=protname_db_reader->getDataByDBKey((char *) id2.c_str());
                        double score = DistanceCalculator::uiLevenshteinDistance( seq1,  seq2);
                        score= 1- (score/std::max(strlen(seq1),strlen(seq2)));
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
            double averagescore;
            if(allagainstall){
                averagescore=(sumofscore / (idswithgo.size()*idswithgo.size()-idswithgo.size()));
            }else{
                averagescore=sumofscore /(idswithgo.size()-1);
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
        if(argc != 7){
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
        bool allagainstall=false;
        if(strcmp(comparisonmode,"yes")==0){
            allagainstall=true;
            Debug(Debug::INFO)<<"all against all comparison";
        }
        Debug(Debug::INFO) << "running keyword evaluation";

        Debug(Debug::INFO) << "Opening clustering database...\n";
        DBReader* cluster_ffindex_reader = new DBReader(cluster_ffindex.c_str(), cluster_ffindex_indexfile.c_str());
        cluster_ffindex_reader->open(DBReader::SORT);
        Debug(Debug::INFO) << "Opening clustering database...\n";
        DBReader* protname_db_reader = new DBReader(keyword_db.c_str(), keyword_indexfile.c_str());
        protname_db_reader->open(DBReader::NOSORT);


        //files
        std::ofstream clusters_full_file;
        clusters_full_file.open(evaluationfolder+fileprefix+"clusters_keywordscores.go");
        clusters_full_file << "algorithm\tclusterid\tid1\tid2\tkeywordscore\n";



        for (int i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


            char *representative=cluster_ffindex_reader->getDbKey(i);
            char *data = cluster_ffindex_reader->getData(i);
            char *idbuffer = new char[255 + 1];
            double sumofscore =0;
            double minscore=1;
            double maxscore=0;

            std::list<std::string> idswithgo;


            //Debug(Debug::INFO) << representative << "\t" << "not available" << "\n";

            while (*data != '\0') {
                Util::parseKey(data, idbuffer);

                if (protname_db_reader->getDataByDBKey(idbuffer) != NULL) {
                    idswithgo.push_back(std::string(idbuffer));


                } else {
                    //Debug(Debug::INFO) << representative << "\t" << idbuffer << "\t" << "not available" <<"\n";
                }

                data = Util::skipLine(data);
            }
            for(std::string id1:idswithgo){
                for(std::string id2:idswithgo){
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
            double averagescore;
            if(allagainstall){
                averagescore=(sumofscore / (idswithgo.size()*idswithgo.size()-idswithgo.size()));
            }else{
                averagescore=sumofscore /(idswithgo.size()-1);
            }


        }




        clusters_full_file.close();

        cluster_ffindex_reader->close();
        cluster_ffindex_reader->~DBReader();
        protname_db_reader->~DBReader();



    }else if (strcmp(argv[1],"-cs")==0) {
        if (argc != 5) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string clusteringfile = argv[2];
        std::string alignmentfile = argv[3];
        std::string outputfile = argv[4];

        convertfiles *cf = new convertfiles();
        cf->getAlignmentscoresForCluster(clusteringfile,alignmentfile,outputfile);

    }else if (strcmp(argv[1],"-df")==0) {
        if (argc != 5) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string domainscorefile = argv[2];
        std::string domainIdentifierFile = argv[3];
        std::string outputfile = argv[4];

        convertfiles *cf = new convertfiles();
        cf->convertDomainFileToFFindex(domainscorefile,domainIdentifierFile,outputfile);

    }else if (strcmp(argv[1],"-ds")==0) {
        if (argc != 7) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string clusteringfile = argv[2];
        std::string alignmentfile = argv[3];
        std::string suffix = argv[4];
        std::string outputfile = argv[5];
        char * comparisonmode=argv[6];
        bool allagainstall=false;
        if(strcmp(comparisonmode,"yes")==0){
            allagainstall=true;
            Debug(Debug::INFO)<<"all against all comparison";
        }
        convertfiles *cf = new convertfiles();
        cf->getDomainScoresForCluster(clusteringfile,alignmentfile,outputfile,suffix,allagainstall);

    }else if (strcmp(argv[1],"-clusterToTsv")==0) {
        if (argc != 5) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string clusteringfile = argv[2];
        std::string suffix = argv[3];
        std::string outputfolder = argv[4];

        convertfiles *cf = new convertfiles();
        cf->convertFfindexToTsv(clusteringfile, suffix, outputfolder);

    }else if(strcmp(argv[1],"-af")==0){
        if (argc != 7) {
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }
        std::string seqDbfile = argv[2];
        std::string alignmentfile = argv[3];
        std::string outputfile = argv[4];
        float threshold = atof(argv[5]);
        int similarityScoreType = atof(argv[6]);

        size_t BUFFER_SIZE = 1000000;
        char* outBuffer = new char[BUFFER_SIZE];
        char *idbuffer1 = new char[255 + 1];
        char *linebuffer = new char[255 + 1];
        char *similarity = new char[255+1];

        DBReader* seqDbr=new DBReader(seqDbfile.c_str(),(seqDbfile+".index").c_str());
        seqDbr->open(DBReader::SORT);
        DBReader* alnDbr=new DBReader(alignmentfile.c_str(),(alignmentfile+".index").c_str());
        seqDbr->open(DBReader::SORT);
        DBWriter* dbw = new DBWriter(outputfile.c_str(), (outputfile+".index").c_str());
        dbw->open();

        for (int i = 0; i < seqDbr->getSize(); i++) {
            char *data = alnDbr->getData(i);
            std::string cluResultsOutString = std::string("");
            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                Util::getLine(data, linebuffer);
                //get similarityscore
                float factor=1;
                float similarityscore;
                if(similarityScoreType==Parameters::APC_ALIGNMENTSCORE){
                    Util::parseByColumnNumber(data, similarity, 1); //column 1 = alignmentscore
                    similarityscore= atof(std::string(similarity).c_str());
                }else if(similarityScoreType==Parameters::APC_COVERAGE){
                    Util::parseByColumnNumber(data, similarity, 2); //column 2 = querycoverage
                    float querycoverage= atof(std::string(similarity).c_str())*factor;
                    Util::parseByColumnNumber(data, similarity, 3); //column 3 = dbcoverage
                    float dbcoverage= atof(std::string(similarity).c_str())*factor;
                    if(querycoverage<dbcoverage){
                        similarityscore=querycoverage;
                    }else{
                        similarityscore=dbcoverage;
                    }

                }else if(similarityScoreType==Parameters::APC_SEQID){
                    Util::parseByColumnNumber(data, similarity, 4); //column 4 = sequence identity
                    similarityscore= atof(std::string(similarity).c_str())*factor;
                }
                else if(similarityScoreType==Parameters::APC_EVAL) {
                    Util::parseByColumnNumber(data, similarity, 5); //column 4 = e value
                    similarityscore = -log(atof(std::string(similarity).c_str()))*factor;
                } else if(similarityScoreType==Parameters::APC_BITSCORE) {
                    Util::parseByColumnNumber(data, similarity, 1); //column 1 = alignmentscore
                    similarityscore= atof(std::string(similarity).c_str());
                    int queryLength=strlen(seqDbr->getDataByDBKey(alnDbr->getDbKey(i)));
                    int dbSeqLength=strlen(seqDbr->getDataByDBKey(idbuffer1));
                    float maxSeqLength=std::max(queryLength,dbSeqLength);

                    //
                    similarityscore= similarityscore/maxSeqLength;

                    //    Debug(Debug::INFO)  << similarityscore <<"\t"<<i<<"\t"<<curr_element<<"\n";
                    //Debug(Debug::INFO)  << similarityscore <<"\n";
                }

                if (similarityscore < threshold) {
                    data = Util::skipLine(data);
                    continue;
                }
                cluResultsOutString=cluResultsOutString+linebuffer+"\n";
                data = Util::skipLine(data);
            }

            // Debug(Debug::INFO)<<data;

            const char* cluResultsOutData = cluResultsOutString.c_str();
            if (BUFFER_SIZE < strlen(cluResultsOutData)){
                Debug(Debug::ERROR) << "Tried to process the alignment list for the query " << alnDbr->getDbKey(i)
                << " , length of the list = " << "\n";
                Debug(Debug::ERROR) << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length()
                << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
                continue;
            }
            memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length()*sizeof(char));
            dbw->write(outBuffer, cluResultsOutString.length(), alnDbr->getDbKey(i));

            //  data = Util::skipLine(data);
            //  }
        }
    }else{
        printHelp();
        Debug(Debug::INFO)<<DistanceCalculator::uiLevenshteinDistance("bla","bla21");
    }




    }

void printHelp() {
    std::string usage("\nEvaluation commands\n");
    usage.append("-go <gofolder> <prot_go_folder> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) >\n");
    usage.append("-pn <prot_name_db> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) >\n");
    usage.append("-kw <keyword_db> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) >\n");
    usage.append("-cs <clustering_file> <alignment_file> <outputfile>\n");
    usage.append("-df <domainscorefile> <domainIdentifierFile> <outputfile>\n");
    usage.append("-ds <clustering_file> <domainscorefile> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) >\n");
    usage.append("-clusterToTsv <clustering_file> <prefix> <outputfolder>\n");
    usage.append("-af  <seqDbfile> <alignmentfile> <outputfile> <cutoff> <scoretype: 1-5>\n");

    Debug(Debug::INFO) << usage << "\n";
    EXIT(EXIT_FAILURE);
}

