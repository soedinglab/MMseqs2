//
// Created by lars on 12.04.15.
//

#include <Util.h>
#include <DistanceCalculator.h>
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
        Debug(Debug::INFO) << "running protein name evaluation";

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



    }

    else{
        Debug(Debug::INFO)<<DistanceCalculator::uiLevenshteinDistance("bla","bla21");
    }




    }

void printHelp() {
    std::string usage("\nEvaluation commands\n");
    usage.append("-go <gofolder> <prot_go_folder> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) >\n");
    usage.append("-pn <prot_name_db> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) >\n");
    usage.append("-kw <keyword_db> <clustering_file> <prefix> <outputfolder> <yes : all against all |no : representative against all(default) >");
    Debug(Debug::INFO) << usage << "\n";
    EXIT(EXIT_FAILURE);
}

