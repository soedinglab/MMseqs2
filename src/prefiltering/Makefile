MAIN_SOURCES := $(shell find ../commons -name "*.cpp")
MAIN_SOURCES += $(shell find . -name "*.cpp")


MAIN_OBJS := $(patsubst %.cpp, %.o, $(MAIN_SOURCES))

CC = g++
CFLAGS = -g -msse4.1  -I../../lib/ffindex/src/ -L../../lib/ffindex/src/ -lffindex  -Wno-write-strings
#CFLAGS = -Wall -Ilib -msse4.1 -m64 -ffast-math -ftree-vectorize -O3 -DOPENMP=1 -fopenmp  -I../../lib/ffindex/src/ -L../../lib/ffindex/src/ -lffindex  -Wno-write-strings 
LDFLAGS =

TARGETS = mmseqs_pref
ifdef HAVE_MPI
CC = mpic++
CFLAGS += -DHAVE_MPI
endif


all: $(TARGETS)

mmseqs_pref: .depend $(MAIN_OBJS)
	$(CC) $(CFLAGS) $(MAIN_OBJS) $(LDFLAGS) -o mmseqs_pref

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
	rm -f .depend *.o
	rm -f $(TARGETS)

.PHONY: clean depend


