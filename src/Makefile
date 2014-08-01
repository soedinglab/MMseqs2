C_FILES := $(shell find commons -name "*.cpp")
C_FILES += $(shell find prefiltering -name "*.cpp" ! -name "Main.cpp")
C_FILES += $(shell find alignment -name "*.cpp" ! -name "Main.cpp")
C_FILES += $(shell find alignment -name "*.C")
C_FILES += $(shell find clustering -name "*.cpp" ! -name "Main.cpp")
C_FILES += workflow/WorkflowFunctions.cpp

PREF_SOURCES := $(C_FILES)
PREF_SOURCES += prefiltering/Main.cpp

ALN_SOURCES := $(C_FILES)
ALN_SOURCES += alignment/Main.cpp

CLU_SOURCES := $(C_FILES)
CLU_SOURCES += clustering/Main.cpp

SEARCH_SOURCES := $(C_FILES)
SEARCH_SOURCES += workflow/Search.cpp

CLUSTER_SOURCES := $(C_FILES)
CLUSTER_SOURCES += workflow/Clustering.cpp

UPDATING_SOURCES := $(C_FILES)
UPDATING_SOURCES += workflow/Updating.cpp 

TT_SOURCES := $(C_FILES)
TT_SOURCES += workflow/TimeTest.h workflow/TimeTest.cpp workflow/TimeTestMain.cpp

PREF_OBJS := $(patsubst %.cpp, %.o, $(PREF_SOURCES))
ALN_OBJS := $(patsubst %.cpp, %.o, $(ALN_SOURCES))
CLU_OBJS := $(patsubst %.cpp, %.o, $(CLU_SOURCES))
SEARCH_OBJS := $(patsubst %.cpp, %.o, $(SEARCH_SOURCES))
CLUSTER_OBJS := $(patsubst %.cpp, %.o, $(CLUSTER_SOURCES))
UPDATING_OBJS := $(patsubst %.cpp, %.o, $(UPDATING_SOURCES))
TT_OBJS := $(patsubst %.cpp, %.o, $(TT_SOURCES))

CC = g++ 
#CFLAGS = -g -Wall  -I../../lib/ffindex/src/ -L../../lib/ffindex/src/ -lffindex  -Wno-write-strings
CFLAGS = -fopenmp -DOPENMP=1 -m64 -ffast-math -ftree-vectorize -O3 -Wno-write-strings -I../lib/ffindex/src/ -L../lib/ffindex/src/ -lffindex -fno-strict-aliasing 
LDFLAGS =

TARGETS = mmseqs_pref mmseqs_aln mmseqs_clu mmseqs_search mmseqs_cluster mmseqs_update time_test

all: $(TARGETS)

mmseqs_pref: .depend $(PREF_OBJS)
	$(CC) $(CFLAGS) $(PREF_OBJS) $(LDFLAGS) -o ../bin/mmseqs_pref

mmseqs_aln: .depend $(ALN_OBJS)
	$(CC) $(CFLAGS) $(ALN_OBJS) $(LDFLAGS) -o ../bin/mmseqs_aln

mmseqs_clu: .depend $(CLU_OBJS)
	$(CC) $(CFLAGS) $(CLU_OBJS) $(LDFLAGS) -o ../bin/mmseqs_clu

mmseqs_search: .depend $(SEARCH_OBJS)
	$(CC) $(CFLAGS) $(SEARCH_OBJS) $(LDFLAGS) -o ../bin/mmseqs_search

mmseqs_cluster: .depend $(CLUSTER_OBJS)
	$(CC) $(CFLAGS) $(CLUSTER_OBJS) $(LDFLAGS) -o ../bin/mmseqs_cluster

mmseqs_update: .depend $(UPDATING_OBJS)
	$(CC) $(CFLAGS) $(UPDATING_OBJS) $(LDFLAGS) -o ../bin/mmseqs_update

time_test: .depend $(TT_OBJS)
	$(CC) $(CFLAGS) $(TT_OBJS) $(LDFLAGS) -o workflow/time_test

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
	rm -f .depend
	rm -f ../bin/mmseqs_pref ../bin/mmseqs_aln ../bin/mmseqs_clu ../bin/mmseqs_search ../bin/mmseqs_cluster ../bin/mmseqs_update workflow/time_test
	rm -f commons/*.o
	rm -f commons/.depend
	rm -f alignment/*.o
	rm -f alignment/.depend
	rm -f prefiltering/*.o
	rm -f prefiltering/.depend
	rm -f clustering/*.o
	rm -f clustering/.depend
	rm -f workflow/*.o
	rm -f workflow/.depend
	rm -f test/*.o
	rm -f test/.depend

.PHONY: all clean depend



