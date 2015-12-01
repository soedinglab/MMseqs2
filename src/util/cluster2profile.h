#ifndef CLUSTER2PROFILE_H
#define CLUSTER2PROFILE_H
#include <string>

int runResult2Profile(std::string resultDb, std::string queryDb, std::string targetDb, std::string outDb,
                      std::string subMatPath, int cpu);

int result2profile(int argn, const char **argv);

#endif
