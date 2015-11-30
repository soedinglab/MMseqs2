//Seg.h & Seg.cpp were modified from the seg package from ncbi ftp.ncbi.nlm.nih.gov:public/seg/seg
#ifndef __SEG_H_
#define __SEG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>

#define LENCORRLIM 120
#define MIN(a,b)        ((a) <= (b) ? (a) : (b))

typedef int Bool;
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

typedef struct alpha0
{
    char *name;
    int size;

    int *index[128];
    char *chars;       /*  [size]  */
    char **charnames;  /*  [size]  */

    Bool caseinvariant;
} Alphabet;

typedef struct seq0
{
    char *seq;
    int start;                     /* for windows */
    int length;

    Bool punctuation;
    Bool seedwin;

    seq0 *parent;       /* for windows */

    int *state;
    int *composition;
    double entropy;

} Sequen;


typedef struct segment0
{
    int begin;
    int end;
    segment0 *next;
} Segment;


class Seg
{
    private:
    int	window;
    int	downset;
    int	upset;
    double	locut;
    double	hicut;
    int	hilenmin;
    int 	overlaps;
    int 	hionly;
    int 	loonly;
    int 	entinfo;
    int 	singleseq;
    int 	prettyseq;
    int 	prettytree;
    int 	charline;
    int 	maxtrim;

    int	thewindow;
    double*	entray;
    char *maskedseq;

    int aaindex[128];
    unsigned char   aaflag[128];
    char aachar[20];

    public:
    Seg(void);
    Seg(int window_inp, int max_seq_len);
    void	initialize(int window_inp);
    ~Seg(void);
    char*	maskseq(const char *seq);
    void	getwin_init(void);
    void	segseq(Sequen *seq, Segment **segs, int offset);
    double*	seqent(Sequen *seq);
    Bool	hasdash(Sequen *win);
    int	findlo(int i, int limit, double *H);
    int	findhi(int i, int limit, double *H);
    void	trim(Sequen *seq, int *leftend, int *rightend);
    double	getprob(int *sv, int total);
    double	lnperm(int *sv, int tot);
    double	lnsass(int *sv);
    void	mergesegs(Sequen *seq, Segment *segs);
    char*	singreport(Sequen *seq, Segment *segs);
    void	appendseg(Segment *segs, Segment *seg);
    void	freesegs(Segment *segs);
    void	closewin(Sequen *win);

    void	entropy_init(int window);
    double	entropy_cal(register int *sv);
    void	closeseq(Sequen *seq);
    Sequen*	openseq(const char *fastaseq);
    Bool	shiftwin1(Sequen *win);
    void	stateon(Sequen *seq);
    void	compon(Sequen *win);
    void	enton(Sequen *win);
    Sequen*	openwin(Sequen *parent, int start, int length);
    int	min(int a, int b);
    void	decrementsv(register int *sv, register int thisclass);
    void	incrementsv(register int *sv, int thisclass);
    double	lnass(int *sv);
    void	upper(register char *string, size_t len);

};

#endif
