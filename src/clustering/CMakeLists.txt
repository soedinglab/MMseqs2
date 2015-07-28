set(clustering_source_files
  AffinityClustering.cpp
  AffinityClustering.h
  Clustering.cpp
  Clustering.h
  LinearMultiArray.h
  Main.cpp
  SetCover.cpp
  SetCover.h
  SetElement.h
  SimpleClustering.cpp
  SimpleClustering.h
  CompareGOTerms.h
  CompareGOTerms.cpp
  DistanceCalculator.h
  SetCover2.cpp
  SetCover2.h
  SetCover3.h
  SetCover3.cpp
  SetCover4.h
  SetCover4.cpp
  AlignmentSymmetry.h
  AlignmentSymmetry.cpp
  SimpleClustering2.h
  SimpleClustering2.cpp
)
add_library(clustering ${clustering_source_files})

target_link_libraries(clustering commons)