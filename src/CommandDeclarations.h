#ifndef COMMANDDECLARATIONS_H
#define COMMANDDECLARATIONS_H

extern int clust(int argc, const char **argv);
extern int search(int argc, const char **argv);
extern int clusteringworkflow(int argc, const char **argv);
extern int clusterupdate(int argc, const char **argv);
extern int prefilter(int argc, const char **argv);
extern int createdb(int argc, const char **argv);
extern int result2flat(int argc, const char **argv);
extern int createindex(int argc, const char **argv);
extern int mergedbs(int argc, const char **argv);
extern int mergeclusters(int argc, const char **argv);
extern int align(int argc, const char **argv);
extern int createseqfiledb(int argc, const char **argv);
extern int swapresults(int argc, const char **argv);
extern int extractorfs(int argc, const char **argv);
extern int convertprofiledb(int argc, const char **argv);
extern int translatenucs(int argc, const char **argv);
extern int result2profile(int argc, const char **argv);
extern int result2msa(int argc, const char **argv);
extern int result2stats(int argc, const char **argv);
extern int splitdb(int argc, const char **argv);
extern int convertalignments(int argc, const char **argv);
extern int createtsv(int argc, const char **argv);
extern int convert2fasta(int argc, const char **argv);
extern int gff2db(int argc, const char **argv);
extern int shellcompletion(int argc, const char** argv);
extern int maskbygff(int argc, const char** argv);
extern int filterdb(int argc, const char** argv);
extern int convertkb(int argc, const char** argv);
extern int subtractdbs(int argc, const char **argv);
extern int computeGOscore(int argc, const char** argv);
extern int clusthash(int argc, const char **argv);
extern int createsubdb(int argc, const char **argv);
extern int summarizeheaders(int argc, const char **argv);
extern int diffseqdbs(int argc, const char **argv);
extern int concatdbs(int argc, const char **argv);
extern int prefixid(int argc, const char** argv);
extern int summarizetabs(int argc, const char **argv);
extern int extractalignedregion(int argc, const char** argv);
extern int extractdomains(int argc, const char **argv);
extern int clustlinear(int argc, const char **argv);

#endif
