add_library(version Version.cpp)

if (VERSION_OVERRIDE)
    target_compile_definitions(version PRIVATE -DGIT_SHA1=${VERSION_OVERRIDE})
else()
    include(GetGitRevisionDescription)
    get_git_head_revision(GIT_REFSPEC GIT_SHA1)
    target_compile_definitions(version PRIVATE -DGIT_SHA1=${GIT_SHA1})
endif()

if (${HAVE_MPI})
    target_compile_definitions(version PRIVATE -DHAVE_MPI=${HAVE_MPI})
endif()
