#ifndef COMMANDDECLARATIONS_H
#define COMMANDDECLARATIONS_H
#include "Command.h"

extern int clust(int argc, const char **argv, const Command& command);
extern int search(int argc, const char **argv, const Command& command);
extern int clusteringworkflow(int argc, const char **argv, const Command& command);
extern int clusterupdate(int argc, const char **argv, const Command& command);
extern int prefilter(int argc, const char **argv, const Command& command);
extern int createdb(int argc, const char **argv, const Command& command);
extern int result2flat(int argc, const char **argv, const Command& command);
extern int createindex(int argc, const char **argv, const Command& command);
extern int mergedbs(int argc, const char **argv, const Command& command);
extern int mergeclusters(int argc, const char **argv, const Command& command);
extern int align(int argc, const char **argv, const Command& command);
extern int createseqfiledb(int argc, const char **argv, const Command& command);
extern int swapresults(int argc, const char **argv, const Command& command);
extern int extractorfs(int argc, const char **argv, const Command& command);
extern int convertprofiledb(int argc, const char **argv, const Command& command);
extern int translatenucs(int argc, const char **argv, const Command& command);
extern int result2profile(int argc, const char **argv, const Command& command);
extern int result2msa(int argc, const char **argv, const Command& command);
extern int result2stats(int argc, const char **argv, const Command& command);
extern int splitdb(int argc, const char **argv, const Command& command);
extern int convertalignments(int argc, const char **argv, const Command& command);
extern int createtsv(int argc, const char **argv, const Command& command);
extern int convert2fasta(int argc, const char **argv, const Command& command);
extern int gff2db(int argc, const char **argv, const Command& command);
extern int shellcompletion(int argc, const char **argv, const Command& command);
extern int maskbygff(int argc, const char **argv, const Command& command);
extern int filterdb(int argc, const char **argv, const Command& command);
extern int convertkb(int argc, const char **argv, const Command& command);
extern int subtractdbs(int argc, const char **argv, const Command& command);
extern int computeGOscore(int argc, const char **argv, const Command& command);
extern int clusthash(int argc, const char **argv, const Command& command);
extern int createsubdb(int argc, const char **argv, const Command& command);
extern int summarizeheaders(int argc, const char **argv, const Command& command);
extern int diffseqdbs(int argc, const char **argv, const Command& command);
extern int concatdbs(int argc, const char **argv, const Command& command);
extern int prefixid(int argc, const char **argv, const Command& command);
extern int summarizetabs(int argc, const char **argv, const Command& command);
extern int extractalignedregion(int argc, const char **argv, const Command& command);
extern int extractdomains(int argc, const char **argv, const Command& command);
extern int clustlinear(int argc, const char **argv, const Command& command);

#endif
