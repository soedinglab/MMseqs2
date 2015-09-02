//
// Created by lars on 12.04.15.
//
#include "CompareGOTerms.h"
#include "Util.h"

#include <dirent.h>
#include <cmath>

CompareGOTerms::CompareGOTerms(std::string go_ffindex,std::string go_ffindex_indexfile,std::string protid_go_ffindex,std::string protid_go_ffindex_indexfile, std::string evaluationfolder) {

    Debug(Debug::INFO) << "Opening GO database...\n";
    go_ffindex_reader = new DBReader(go_ffindex.c_str(), go_ffindex_indexfile.c_str());
    go_ffindex_reader->open(DBReader::NOSORT);
    Debug(Debug::INFO) << "Opening Protein to GO database...\n";
    protid_go_ffindex_reader = new DBReader(protid_go_ffindex.c_str(), protid_go_ffindex_indexfile.c_str());
    protid_go_ffindex_reader->open(DBReader::NOSORT);
    if (evaluationfolder.back() != '/') {
        evaluationfolder = evaluationfolder + '/';
    }
    DIR *dir = opendir(evaluationfolder.c_str());
    if (dir) {
        closedir(dir);
    } else {
        Debug(Debug::ERROR) << "Evaluationfolder: "<<evaluationfolder<<" does not exist!\n";
        EXIT(EXIT_FAILURE);
    }
    this->evaluationfolder = evaluationfolder;

}

CompareGOTerms::~CompareGOTerms() {
    go_ffindex_reader->~DBReader();
    protid_go_ffindex_reader->~DBReader();
    delete[] is_a_relation;
    delete[] is_a_relation_data;
    delete[] is_a_relation_size;
    //counts of goterms
    delete[] count_goterm;

    delete[]count_accumulated_goterm;

    delete[] parentsets;

//    std::map<int, int> goterm_to_index;
  //  std::map<int, int> index_togoterm;
}

void CompareGOTerms::init() {


    //buffers
    char *singletermbuffer = new char[11];
    char *parentterm = new char[255 + 1];
//constants
    total_go_number = go_ffindex_reader->getSize();
    size_t m = protid_go_ffindex_reader->getSize();

    //calculate number of pairs for initialisation of array
    is_a_relation_number = (go_ffindex_reader->getDataSize() - total_go_number * 2) / 10;
    is_a_relation_data = new int[is_a_relation_number];
    is_a_relation = new int *[total_go_number];
    is_a_relation_size = new int[total_go_number];

    count_goterm = new int[total_go_number];
    count_accumulated_goterm = new int[total_go_number];
    //first iteration, create mapping to index
    for (size_t i = 0; i < total_go_number; i++) {
        //initialisierung
        count_goterm[i]=0;
        count_accumulated_goterm[i]=0;
        // char *data = go_ffindex_reader->getData(i);
        char *childterm = go_ffindex_reader->getDbKey(i);
        strncpy(singletermbuffer, childterm + 3, 7);
        singletermbuffer[7] = '\0';
        index_togoterm[i] = atof(std::string(singletermbuffer).c_str());
        goterm_to_index[atof(std::string(singletermbuffer).c_str())] = i;
        // Debug(Debug::INFO) << childterm << "\t" << i << "\t" << index_togoterm[i] << "\t" <<
        // convert_index_toGOterm(i) << "\t" << convert_GOterm_to_index(convert_index_toGOterm(i)) << "\n";

    }

    //second iterationgetGOListforProtein
    size_t is_a_relation_current_position = 0;
    for (size_t i = 0; i < total_go_number; i++) {
        char *data = go_ffindex_reader->getData(i);
        is_a_relation[i] = is_a_relation_data + is_a_relation_current_position;
        is_a_relation_size[i] = 0;
        while (*data != '\0') {
            Util::parseKey(data, parentterm);
            int counter = 0;
            singletermbuffer = new char[11];
            char *position = parentterm;
            while (parentterm[counter] != '\0') {
                strncpy(singletermbuffer, position + 3, 7);
                singletermbuffer[10] = '\0';
                //  Debug(Debug::INFO) << singletermbuffer << "\t";
                counter += 10;
                position = position + 10;
                if(atof(std::string(singletermbuffer).c_str())==0){

                }else{
                    is_a_relation_data[is_a_relation_current_position++] = convert_GOterm_to_index(
                            atof(std::string(singletermbuffer).c_str()));
                    is_a_relation_size[i]++;
                }

            }
            //Debug(Debug::INFO) << "\n";
            //


            data = Util::skipLine(data);
        }
    }
    for (size_t i = 0; i < total_go_number; i++) {
        //Debug(Debug::INFO) << convert_index_toGOterm(i)<<":\t";
        for (size_t j = 0; j < is_a_relation_size[i]; j++) {
            //   Debug(Debug::INFO) << is_a_relation[i][j]<<"\t";
        }
    }
    //bestimme vaterknoten

    //Debug(Debug::INFO) << <<"\t";
    parentsets = new std::set<int>[total_go_number];
    for (size_t i = 0; i < total_go_number; i++) {
        parentsets[i] = *new std::set<int>();
        compute_parentnodes(parentsets[i], i);
    }
    /*    //print
   for (size_t i = 0; i < n; i++) {
       Debug(Debug::INFO) << convert_index_toGOterm(i)<<":\n";
       for(int parentid : parentsets[i]){
           Debug(Debug::INFO) << convert_index_toGOterm(parentid)<<"\t";
       }
       Debug(Debug::INFO) << "\n";

   }*/

    //count occurence in the dataset

    for (size_t j = 0; j < m; j++) {

        char *data = protid_go_ffindex_reader->getData(j);
        char *protid = protid_go_ffindex_reader->getDbKey(j);
        size_t counter = 0;
        count_goterm[convert_GOterm_to_index(atof(strncpy(singletermbuffer, data + counter + 3, 7)))]++;
        while (data[counter] != '\0') {
            if (data[counter] == '.' && (data[counter + 1] != '\0' && data[counter + 1] != '\n')) {
                count_goterm[convert_GOterm_to_index(atof(strncpy(singletermbuffer, data + counter + 4, 7)))]++;
            }
            counter++;
        }


    }
    //transfer counts to parents
    count_goterm_total_sum=0.0;
    for (size_t i = 0; i < total_go_number; i++) {
        count_goterm_total_sum+= count_goterm[i];
        // Debug(Debug::INFO) << convert_index_toGOterm(i) << ":\n";
        for (int parentid : parentsets[i]) {
            //Debug(Debug::INFO) << convert_index_toGOterm(parentid) <<"\t";
            count_accumulated_goterm[parentid] += count_goterm[i];
        }
       // Debug(Debug::INFO) << "\n";

    }
    //print

    for (size_t i = 0; i < total_go_number; i++) {
        //Debug(Debug::INFO) << convert_index_toGOterm(i) << "\t" << count_goterm[i] << "\t" << count_accumulated_goterm[i] <<"\n";
    }
    //Debug(Debug::INFO) <<count_goterm_total_sum;

delete [] singletermbuffer;
delete[] parentterm;
}

double CompareGOTerms::compare_protein_ids(char *prot1, char *prot2) {
    return similarity_of_list(getGOListforProtein(prot1),getGOListforProtein(prot2));
}

void CompareGOTerms::all_against_all_comparison() {
    for (size_t i = 0; i < total_go_number; i++) {
        Debug(Debug::INFO) <<convert_index_toGOterm(i) << "\t";
        for (int j = 0; j < total_go_number; ++j) {
            if(count_goterm[i]>0&&count_goterm[j]>0) {
                Debug(Debug::INFO) << similarity(i, j) << "\t";
            }
        }
        Debug(Debug::INFO) <<"\n";
    }
}
void CompareGOTerms::all_against_all_comparison_proteinset() {
    for (size_t i = 0; i < count_goterm_total_sum; i++) {
        char*protein1=protid_go_ffindex_reader->getDbKey(i);
        Debug(Debug::INFO) <<protein1 << "\t";
        for (int j = 0; j < total_go_number; ++j) {
            char*protein2=protid_go_ffindex_reader->getDbKey(j);

                Debug(Debug::INFO) << compare_protein_ids(protein1, protein2) << "\t";

        }
        Debug(Debug::INFO) <<"\n";
    }
}


double CompareGOTerms::similarity_of_list(std::list<int> ids1, std::list<int> ids2) {
    double result = 0;
    for (int id1:ids1) {
        double maxscore = 0;
        for (int id2:ids2) {
            maxscore = std::max(maxscore, similarity(id1, id2));
        }
        result += maxscore;
    }
    for (int id2:ids2) {
        double maxscore = 0;
        for (int id1:ids1) {
            maxscore = std::max(maxscore, similarity(id1, id2));
        }
        result += maxscore;
    }
    return result / (ids1.size() + ids2.size());
}

double CompareGOTerms::similarity(int id1, int id2) {
            return (2*std::log(most_specific_parent(id1,id2)/count_goterm_total_sum))/(std::log(count_accumulated_goterm[id1]/count_goterm_total_sum)+std::log(count_accumulated_goterm[id2]/count_goterm_total_sum));
}

int CompareGOTerms::most_specific_parent(int id1, int id2) {
    int result = count_goterm_total_sum;
    for (int parentid1 : parentsets[id1]) {
        for (int parentid2 : parentsets[id2]) {
            if(parentid1==parentid2){
                result=std::min(result,count_accumulated_goterm[parentid1]);
            }
        }
    }
    return result;
}

std::list<int> CompareGOTerms::getGOListforProtein(char* protid) {
    char *singletermbuffer = new char[11];
    char *data = protid_go_ffindex_reader->getDataByDBKey(protid);
    std::list<int> result;
    if(NULL == data){ // check if file contains entry
        delete [] singletermbuffer;
        return result;
    }
    size_t counter = 0;
    result.push_back(convert_GOterm_to_index(atof(strncpy(singletermbuffer, data + counter + 3, 7))));
    while (data[counter] != '\0') {
        if (data[counter] == '.' && (data[counter + 1] != '\0' && data[counter + 1] != '\n')) {
            result.push_back(convert_GOterm_to_index(atof(strncpy(singletermbuffer, data + counter + 4, 7))));
        }
        counter++;
    }
    delete [] singletermbuffer;
    return result;
}

int CompareGOTerms::convert_GOterm_to_index(int Goterm){
    return goterm_to_index[Goterm];
}
int CompareGOTerms::convert_index_toGOterm(int index){
    return index_togoterm[index];
}

void  CompareGOTerms::compute_parentnodes(std::set<int>& result,int id){
            result.insert(id);
    if(is_a_relation_size[id]==0){
        return;
    }else{
        for (int i = 0; i < is_a_relation_size[id]; ++i) {
            compute_parentnodes(result,is_a_relation[id][i]);
        }
        return;
    }
    
    
}


void CompareGOTerms::run_evaluation_mmseqsclustering(std::string cluster_ffindex,
                                                     std::string cluster_ffindex_indexfile,
                                                     std::string fileprefix,
                                                     std::string filesuffix, bool allagainstall) {
    Debug(Debug::INFO) << "Opening clustering database...\n";
    DBReader* cluster_ffindex_reader = new DBReader(cluster_ffindex.c_str(), cluster_ffindex_indexfile.c_str());
    cluster_ffindex_reader->open(DBReader::SORT);
    //files
    std::ofstream clusters_summary_file;
    clusters_summary_file.open(evaluationfolder+fileprefix+"clusters_summary.go"+filesuffix);
    clusters_summary_file << "algorithm\tgocategory\tclusterid\tclustersize\twithgo\tmissinggo\tavggoscore\tmingoscore\tmaxgoscore\n";

    std::ofstream clusters_full_file;
    clusters_full_file.open(evaluationfolder+fileprefix+"clusters_allscores.go"+filesuffix);
    clusters_full_file << "algorithm\tgocategory\tclusterid\tid1\tid2\tgoscore\n";

    std::ofstream summary_file;
    summary_file.open(evaluationfolder+fileprefix+"summary.go"+filesuffix);
    summary_file <<"algorithm\tgocategory\tclusternumber\tnwithgo\tnmissinggo\toccurenceofgoterms\n";

    std::ofstream binned_scores_file;
    binned_scores_file.open(evaluationfolder+fileprefix+"binned_scores.go"+filesuffix);
    binned_scores_file <<"algorithm\tgocategory\tgoscore\tgoavg\tgomin\tgomax\n";
    int clusterwithgo=0;
    int clusterwithoutgo=0;
    int binsize=10;
    int *binned_avg=new int[binsize];
    int *binned_min=new int[binsize];
    int *binned_max=new int[binsize];
    for (int j = 0; j < binsize; ++j) {
        binned_avg[j]=0;
        binned_min[j]=0;
        binned_max[j]=0;
    }
    for (int i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


        char *representative=cluster_ffindex_reader->getDbKey(i);
        char *data = cluster_ffindex_reader->getData(i);
        char *idbuffer = new char[255 + 1];
        int withgo=0;
        int withoutgo=0;
        double sumofscore =0;
        double minscore=1;
        double maxscore=0;

        std::list<std::string> idswithgo;


            //Debug(Debug::INFO) << representative << "\t" << "not available" << "\n";

            while (*data != '\0') {
                Util::parseKey(data, idbuffer);

                    if (protid_go_ffindex_reader->getDataByDBKey(idbuffer) != NULL) {
                        idswithgo.push_back(std::string(idbuffer));

                        withgo++;
               clusterwithgo=0;
            clusterwithoutgo=0;
                    } else {
                        //Debug(Debug::INFO) << representative << "\t" << idbuffer << "\t" << "not available" <<"\n";
                        withoutgo++;
                    }

                data = Util::skipLine(data);
            }
            for(std::string id1:idswithgo){
                for(std::string id2:idswithgo){
                    if (std::string(id1) != std::string(id2)) {
                        double score = compare_protein_ids((char *) id1.c_str(), (char *) id2.c_str());
                        sumofscore += score;
                        minscore=std::min(score,minscore);
                        maxscore=std::max(score,maxscore);
                        clusters_full_file << fileprefix << "\t"<< filesuffix<< "\t"<< representative << "\t" << id1 << "\t" << id2 << "\t" << score << "\n";
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
            if(idswithgo.size()>1){
                binned_avg[(int)(averagescore*binsize)%binsize]++;
                binned_min[(int)(minscore*binsize)%binsize]++;
                binned_max[(int)(std::min(maxscore,0.999)*binsize)%binsize]++;
                clusterwithgo++;
            }else{
                clusterwithoutgo++;
            }
       // if(idswithgo.size()>0) {
            clusters_summary_file << fileprefix << "\t" << filesuffix << "\t" << representative << "\t" <<
                                                   withgo + withoutgo << "\t" << withgo << "\t" << withoutgo << "\t" <<
                    averagescore << "\t" << minscore << "\t" << maxscore << "\n";
        //}
    }
    for (int j = 0; j < binsize; ++j) {
        binned_scores_file << fileprefix << "\t"<< filesuffix<< "\t"<< (j/(double)binsize)<<"\t"<<binned_avg[j]<<"\t"<<binned_min[j]<<"\t"<<binned_max[j]<<"\n";
    }

    summary_file << fileprefix << "\t"<< filesuffix<< "\t"<< cluster_ffindex_reader->getSize()<<"\t"<<clusterwithgo <<"\t"<<clusterwithoutgo <<"\t"<<this->count_goterm_total_sum <<"\n";

    clusters_summary_file.close();
    clusters_full_file.close();
    summary_file.close();
    binned_scores_file.close();
    cluster_ffindex_reader->close();
    cluster_ffindex_reader->~DBReader();

}

