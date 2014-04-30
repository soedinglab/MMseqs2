MAIN_SOURCES := $(shell find ../src/commons -name "*.cpp")
MAIN_SOURCES += $(shell find ../src/prefiltering -name "*.cpp" ! -name "Main.*")


MAIN_OBJS := $(patsubst %.cpp, %.o, $(MAIN_SOURCES))

CC = g++
CFLAGS = -ggdb -msse4.1  -I../lib/ffindex/src/ -L../lib/ffindex/src/ -lffindex  -Wno-write-strings
#CFLAGS = -g -Wall -msse4.1 -m64 -ffast-math -ftree-vectorize -fstrict-aliasing -O3 -DOPENMP=1 -fopenmp  -I../lib/ffindex/src/ -L../lib/ffindex/src/ -lffindex  -Wno-write-strings 
CFLAGS += -I../src/commons -I../src/prefiltering 
TARGETS = ffindex_database_merge 
TARGETS += ffindex_from_fasta
TARGETS += mmseqs_pref_index
TARGETS += ffindex_to_fasta
ifdef HAVE_MPI
CC = mpic++
CFLAGS += -DHAVE_MPI
endif


all: $(TARGETS)
ffindex_database_merge: .depend ffindex_database_merge.o $(MAIN_OBJS)
	$(CC) $(CFLAGS) $(MAIN_OBJS) ffindex_database_merge.o -o ffindex_database_merge

ffindex_database_merge.o: ffindex_database_merge.cpp
	$(CC) $(CFLAGS) -c ffindex_database_merge.cpp -o ffindex_database_merge.o

ffindex_to_fasta: .depend ffindex_to_fasta.o $(MAIN_OBJS)
	$(CC) $(CFLAGS) $(MAIN_OBJS) ffindex_to_fasta.o -o ffindex_to_fasta

ffindex_to_fasta.o: ffindex_to_fasta.cpp
	$(CC) $(CFLAGS) -c ffindex_to_fasta.cpp -o ffindex_to_fasta.o

mmseqs_pref_index: .depend mmseqs_pref_index.o $(MAIN_OBJS)
	$(CC) $(CFLAGS) $(MAIN_OBJS) mmseqs_pref_index.o -o mmseqs_pref_index
mmseqs_pref_index.o: mmseqs_pref_index.o
	$(CC) $(CFLAGS) -c mmseqs_pref_index.cpp -o mmseqs_pref_index.o

ffindex_from_fasta: ffindex_from_fasta .depend 
	gcc $(CFLAGS) -std=c99 ffindex_from_fasta.c -o ffindex_from_fasta

depend: .depend

.depend: cmd = gcc -MM -MF depend $(var); cat depend >> .depend;
.depend:
	@echo "Generating dependencies..."
	@$(foreach var, $(C_FILES), $(cmd))
	@rm -f depend

-include .depend

# These are the pattern matching rules. In addition to the automatic
# variables used here, the variable $* that matches whatever % stands for
# can be useful in special cases.
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%: %.cpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm ../src/commons/*.o
	rm ../src/prefiltering/*.o
	rm -f .depend *.o
	rm -f $(TARGETS)

.PHONY: clean depend


