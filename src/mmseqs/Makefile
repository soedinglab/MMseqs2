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
#CFLAGS = -ggdb -msse4.1 $(INCLUDES) -L../../lib/ffindex/src/ -lffindex  -Wno-write-strings -lz
CFLAGS = -Wall -msse4.1 -m64 -ffast-math -ftree-vectorize -fstrict-aliasing -O3 -DOPENMP=1 -fopenmp $(INCLUDES) -L../../lib/ffindex/src/ -lffindex -lz  -Wno-write-strings 
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
	rm ../commons/*.o
	rm ../prefiltering/*.o
	rm ../alignment/*.o
	rm ../clustering/*.o
	rm ../workflow/*.o
	rm ../util/*.o
	rm -f *.o
	rm -f $(TARGETS)



