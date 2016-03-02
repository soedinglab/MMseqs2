set(clustering_source_files
  Clustering.cpp
  Clustering.h
  Main.cpp
  SetElement.h
  CompareGOTerms.h
  CompareGOTerms.cpp
  DistanceCalculator.h
  ClusteringAlgorithms.h
  ClusteringAlgorithms.cpp
  AlignmentSymmetry.h
  AlignmentSymmetry.cpp
)
add_library(clustering ${clustering_source_files})

target_link_libraries(clustering commons)
