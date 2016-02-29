set(workflow_source_files
  Clustering.cpp
  Search.cpp
  Updating.cpp
)
add_library(workflow ${workflow_source_files})

target_link_libraries(workflow prefiltering alignment util clustering)
