include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

add_library(version Version.cpp)
target_compile_definitions(version PRIVATE -DGIT_SHA1=${GIT_SHA1})
