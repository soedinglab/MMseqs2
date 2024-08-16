#ifndef CLUSTER_SPECIES_H
#define CLUSTER_SPECIES_H

#include <vector>
#include <unordered_map>


class ClusterSpecies {
public:
    
    // std::unordered_map<unsigned int, std::unordered_map<unsigned int, std::vector<unsigned int>>> clusterEntriesPerSpecies;
    // // std::unordered_map<unsigned int, size_t> speciesClusterMemCounts;
    // void addProtein(unsigned int speciesId, unsigned int clusterRepId, unsigned int proteinId);
    // void removeSpecies(unsigned int speciesId);
    // int findRepSpecies();

    // template <typename T>
    // struct ProteomeLength{
    //     T id;
    //     unsigned int length;
    // };
    
    // std::vector<ProteomeLength<unsigned int>> proteomeAAToatalLength;

    // void addProteomeLength(unsigned int id, unsigned int length){
    //     ProteomeLength<unsigned int> pl;
    //     pl.id = id;
    //     pl.length = length;
    //     proteomeAAToatalLength.push_back(pl);
    // }

    // unsigned int getProteomeLength(unsigned int id){
    //     for (auto const& pl : proteomeAAToatalLength){
    //         if (pl.id == id){
    //             return pl.length;
    //         }
    //     }
    //     return 0;
    // }
};

#endif //CLUSTER_SPECIES_H