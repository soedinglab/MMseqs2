set(workflow_source_files
        Cluster.cpp
  Search.cpp
        ClusterUpdate.cpp
)
add_library(workflow ${workflow_source_files})

target_link_libraries(workflow prefiltering alignment util clustering)
