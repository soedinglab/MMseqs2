MAIN_SOURCES := $(shell find . -name "*.cpp")
MAIN_SOURCES += $(shell find ../alignment -name "*.cpp")
MAIN_SOURCES += $(shell find ../commons -name "*.cpp")
MAIN_SOURCES += $(shell find ../prefiltering -name "*.cpp")
MAIN_SOURCES += $(shell find ../clustering -name "*.cpp")
MAIN_SOURCES += $(shell find ../workflow -name "*.cpp")
MAIN_SOURCES += $(shell find ../util -name "*.cpp")

HEADER_SOURCES := $(shell find . -name "*.h")
HEADER_SOURCES += $(shell find ../alignment -name "*.h")
HEADER_SOURCES += $(shell find ../commons -name "*.h")
HEADER_SOURCES += $(shell find ../prefiltering -name "*.h")
HEADER_SOURCES += $(shell find ../clustering -name "*.h")
HEADER_SOURCES += $(shell find ../workflow -name "*.h")
HEADER_SOURCES += $(shell find ../util -name "*.h")

MAIN_OBJS := $(patsubst %.cpp, %.o, $(MAIN_SOURCES))


INCLUDES=-I../commons/ -I../../lib/ffindex/src/ -I../../lib/kseq/ -I../prefiltering/ -I../alignment/ -I../clustering/ -I../workflow/ -I../util/
CC = g++
#CFLAGS = -std=c++11 -ggdb -msse4.2 $(INCLUDES) -L../../lib/ffindex/src/ -lffindex  -Wno-write-strings -lz
CFLAGS = -std=c++11 -g -Wall -msse4.2 -m64 -ffast-math -ftree-vectorize -fno-strict-aliasing -O3 -DOPENMP=1 -fopenmp $(INCLUDES) -L../../lib/ffindex/src/ -lffindex -lz  -Wno-write-strings 
LDFLAGS =

TARGETS = mmseqs
ifdef HAVE_MPI
CC = mpic++
CFLAGS += -DHAVE_MPI
endif


all: $(TARGETS)

mmseqs: $(MAIN_OBJS) 
	$(CC) $(CFLAGS) $(MAIN_OBJS) $(LDFLAGS) -o mmseqs

# These are the pattern matching rules. In addition to the automatic
# variables used here, the variable $* that matches whatever % stands for
# can be useful in special cases.
%.o: %.cpp ${HEADER_SOURCES}
	$(CC) $(CFLAGS) -c $< -o $@

%: %.cpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f ../commons/*.o
	rm -f ../prefiltering/*.o
	rm -f ../alignment/*.o
	rm -f ../clustering/*.o
	rm -f ../workflow/*.o
	rm -f ../util/*.o
	rm -f *.o
	rm -f $(TARGETS)



