C_FILES := $(shell find ../commons -name "*.cpp")
C_FILES += $(shell find ../prefiltering -name "*.cpp" ! -name "Main.cpp")
C_FILES += $(shell find ../alignment -name "*.cpp" ! -name "Main.cpp")
C_FILES += $(shell find ../alignment -name "*.C")
C_FILES += $(shell find ../clustering -name "*.cpp" ! -name "Main.cpp")

CC_SOURCES := $(C_FILES)
CC_SOURCES += CascadedClustering.cpp

UPDATING_SOURCES := $(C_FILES)
UPDATING_SOURCES += Updating.cpp 

TT_SOURCES := $(C_FILES)
TT_SOURCES += TimeTest.h TimeTest.cpp TimeTestMain.cpp

CC_OBJS := $(patsubst %.cpp, %.o, $(CC_SOURCES))
UPDATING_OBJS := $(patsubst %.cpp, %.o, $(UPDATING_SOURCES))
TT_OBJS := $(patsubst %.cpp, %.o, $(TT_SOURCES))

CC = g++ 
#CFLAGS = -g  -I../../lib/ffindex/src/ -L../../lib/ffindex/src/ -lffindex  -Wno-write-strings
CFLAGS = -Wall -fopenmp -DOPENMP=1 -m64 -ffast-math -ftree-vectorize -O3 -Wno-write-strings -I../../lib/ffindex/src/ -L../../lib/ffindex/src/ -lffindex -fno-strict-aliasing 
LDFLAGS =

TARGETS = cascaded_clustering updating time_test

all: $(TARGETS)

cascaded_clustering: .depend $(CC_OBJS)
	$(CC) $(CFLAGS) $(CC_OBJS) $(LDFLAGS) -o cascaded_clustering

updating: .depend $(UPDATING_OBJS)
	$(CC) $(CFLAGS) $(UPDATING_OBJS) $(LDFLAGS) -o updating

time_test: .depend $(TT_OBJS)
	$(CC) $(CFLAGS) $(TT_OBJS) $(LDFLAGS) -o time_test

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

.PHONY: all clean depend



