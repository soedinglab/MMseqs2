#ifndef COMMANDDECLARATIONS_H
#define COMMANDDECLARATIONS_H

extern int cluster(int argc, const char **argv);
extern int search(int argc, const char **argv);
extern int clusteringworkflow(int argc, const char **argv);
extern int clusterupdate(int argc, const char **argv);
extern int prefilter(int argc, const char **argv);
extern int createdb(int argc, const char **argv);
extern int createfasta(int argc, const char **argv);
extern int createindex(int argc, const char **argv);
extern int mergeffindex(int argc, const char **argv);
extern int mergecluster(int argc, const char **argv);
extern int alignment(int argc, const char **argv);
extern int addsequences(int argc, const char **argv);
extern int swapresults(int argc, const char **argv);
extern int extractorf(int argc, const char **argv);
extern int createprofiledb(int argc, const char **argv);
extern int translatenucleotide(int argc, const char **argv);
extern int timetest(int argc, const char **argv);
extern int result2profile(int argc, const char **argv);
extern int result2msa(int argc, const char **argv);
extern int splitffindex(int argc, const char **argv);
extern int formatalignment(int argc, const char **argv);
extern int createtsv(int argc, const char **argv);
extern int rebuildfasta(int argc, const char **argv);
extern int gff2ffindex(int argc, const char **argv);
extern int shellcompletion(int argc, const char** argv);
extern int maskbygff(int argc, const char** argv);
extern int filterdb(int argc, const char** argv);
extern int convertkb(int argc, const char** argv);
extern int substractresult(int argc, const char** argv);
extern int result2newick(int argc, const char** argv);
extern int kbtotsv(int argc, const char** argv);
extern int computeGOscore(int argc, const char** argv);
extern int detectredundancy(int argc, const char** argv);
extern int order(int argc, const char** argv);
extern int summarizeheaders(int argc, const char **argv);
extern int diff(int argc, const char** argv);
extern int dbconcat(int argc, const char** argv);
extern int prefixid(int argc, const char** argv);
extern int summarizetabs(int argc, const char **argv);
extern int extractalignedregion(int argc, const char** argv);
extern int extractdomains(int argc, const char **argv);

#endif
