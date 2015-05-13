set(workflow_source_files
  Clustering.cpp
  Search.cpp
  TimeTest.cpp
  TimeTestMain.cpp
  Updating.cpp
  WorkflowFunctions.cpp
  WorkflowFunctions.h
)
add_library(workflow ${workflow_source_files})

target_link_libraries(workflow prefiltering alignment util clustering)
