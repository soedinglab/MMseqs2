MAIN_SOURCES := $(shell find ../commons -name "*.cpp")
MAIN_SOURCES += $(shell find ../prefiltering -name "*.cpp" ! -name "Main.cpp")
MAIN_SOURCES += $(shell find ../alignment -name "*.cpp" ! -name "Main.cpp")
TARGETS := $(shell find . -name "*.cpp")
TARGETS := $(patsubst %.cpp, %, $(TARGETS))

MAIN_OBJS := $(patsubst %.cpp, %.o, $(MAIN_SOURCES))

CC = g++
#CFLAGS = -g -pg  -I../../lib/ffindex/src/ -I../commons/ -I../prefiltering/ -L../../lib/ffindex/src/ -lffindex  -Wno-write-strings
CFLAGS = -Wall -std=c++11 -Ilib -m64 -ffast-math -ftree-vectorize -fstrict-aliasing -O3 -DOPENMP=1 -fopenmp -I../commons/ -I../prefiltering/  -I../../lib/ffindex/src/ -L../../lib/ffindex/src/ -Wno-write-strings 
LDFLAGS =

all: $(TARGETS)

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
	$(CC) $(CFLAGS) $(MAIN_OBJS) ../../lib/ffindex/src/libffindex.a $(LDFLAGS) -o $@ $<

clean:
	rm -f .depend *.o
	rm -f $(TARGETS)

.PHONY: clean depend



