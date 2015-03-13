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
)
add_library(clustering ${clustering_source_files})

target_link_libraries(clustering commons)