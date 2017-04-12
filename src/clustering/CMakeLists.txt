set(clustering_source_files
  Clustering.cpp
  Clustering.h
  Main.cpp
  SetElement.h
  DistanceCalculator.h
  DistanceCalculator.cpp
  ClusteringAlgorithms.h
  ClusteringAlgorithms.cpp
  AlignmentSymmetry.h
  AlignmentSymmetry.cpp
)
add_library(clustering ${clustering_source_files})

target_link_libraries(clustering commons)
