#include "ClusterSpecies.h"
void ClusterSpecies::addProtein(unsigned int speciesId, unsigned int clusterRepId, unsigned int proteinId) {
   clusterEntriesPerSpecies[speciesId][clusterRepId].push_back(proteinId);
//    speciesClusterMemCounts[speciesId]++;
}


void ClusterSpecies::removeSpecies(unsigned int speciesId) {
    clusterEntriesPerSpecies.erase(speciesId);
    // speciesClusterMemCounts.erase(speciesId);
}

int ClusterSpecies::findRepSpecies(){
    //find max entry in speciesClusterMemCounts
    unsigned int maxSpeciesId = 0;
    size_t maxSpeciesSize = 0;
    for (auto const& species : clusterEntriesPerSpecies){
        size_t speiciesSize = species.second.size();
        if (speiciesSize > maxSpeciesSize){
            maxSpeciesId = species.first;
        }
    }

    return maxSpeciesId;
    // for (auto const& x : speciesClusterMemCounts)
    // {
    //     if (x.second > maxSpeciesCount){
    //         maxSpeciesCount = x.second;
    //         maxSpeciesId = x.first;
    //     }
    // }
}