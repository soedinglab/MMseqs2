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
extern int alignment(int argc, const char **argv);
extern int clusteringtofastadb(int argc, const char **argv);
extern int swapresults(int argc, const char **argv);
extern int extractorf(int argc, const char **argv);
extern int createprofiledb(int argc, const char **argv);
extern int translatenucleotide(int argc, const char **argv);
extern int timetest(int argc, const char **argv);
extern int legacycs219(int argc, const char **argv);

#endif
