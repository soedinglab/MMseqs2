INSTALL_DIR?=..
libdir=`([ -d /usr/lib64 ] && echo lib64) || echo lib`
LANG=C
CC=gcc -fopenmp
#CC=icc -openmp
CFLAGS+=-std=c99 -O2 -g -fPIC -Wall -pedantic
LDFLAGS=-lffindex -L.
ifdef FFINDEX_STATIC
	LDFLAGS+=-static
	FFINDEX_OBJECTS=libffindex.a
else
	FFINDEX_OBJECTS=libffindex.so libffindex.so.0.1 
endif

#-Wl,--version-script=libffindex.map

OBJECTS=ffindex_apply ffindex_from_fasta ffindex_get ffindex_build ffindex_modify libffindex.a libffindex.so ffindex_unpack

ifdef HAVE_MPI
OBJECTS+=ffindex_apply_mpi
endif

all: $(OBJECTS)

ffindex.o: ffindex.c ffindex.h ffutil.h

ffindex_apply: $(FFINDEX_OBJECTS) ffindex_apply.o ffindex.h

ffindex_apply_mpi: $(FFINDEX_OBJECTS) ffindex_apply_mpi.o ffindex.h
	mpicc $(CFLAGS) -o $@  ffindex_apply_mpi.o $(LDFLAGS)

ffindex_apply_mpi.o: ffindex_apply_mpi.c ffindex.h
	mpicc $(CFLAGS) -c ffindex_apply_mpi.c -o $@ 

ffindex_from_fasta: $(FFINDEX_OBJECTS) ffindex_from_fasta.o ffindex.h
	$(CC) -o $@ ffindex_from_fasta.o $(LDFLAGS) -L. -lffindex

ffindex_get: $(FFINDEX_OBJECTS) ffindex_get.o ffindex.h
	$(CC) -o $@ ffindex_get.o $(LDFLAGS) -L. -lffindex

ffindex_build: $(FFINDEX_OBJECTS) ffindex_build.o ffindex.h
	$(CC) -o $@ ffindex_build.o $(LDFLAGS) -L. -lffindex

ffindex_modify: $(FFINDEX_OBJECTS) ffindex_modify.o ffindex.h
	$(CC) -o $@ ffindex_modify.o $(LDFLAGS) -L. -lffindex

ffindex_unpack: $(FFINDEX_OBJECTS) ffindex_unpack.o ffindex.h
	$(CC) -o $@ ffindex_unpack.o $(LDFLAGS) -L. -lffindex

libffindex.so.0.1: ffindex.o ffutil.o
	$(CC) -Wl,-z,defs -shared -Wl,-soname,libffindex.so.0.1 -o libffindex.so.0.1 ffindex.o ffutil.o -lc

libffindex.so: libffindex.so.0.1
	ln -sf libffindex.so.0.1 libffindex.so

libffindex.a: ffindex.o ffutil.o
	ar rcs libffindex.a ffutil.o ffindex.o 

test:
	rm -f /tmp/test.data /tmp/test.ffindex
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" ./ffindex_build -s /tmp/test.data /tmp/test.ffindex ../test/data ../test/data2
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" ./ffindex_get  /tmp/test.data /tmp/test.ffindex a b foo | tee /tmp/test.out
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" ./ffindex_modify -u /tmp/test.ffindex b
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" ./ffindex_get  /tmp/test.data /tmp/test.ffindex a b foo | tee /tmp/test-unlink.out
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" ./ffindex_build -a -s /tmp/test.data /tmp/test.ffindex ../test/testfile
	@echo ======== compare results =========
	diff /tmp/test.out test.should
	diff /tmp/test-unlink.out test-unlink.should
	@echo ======== error handling ==========
	@echo
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" ./ffindex_get  /tmp/test.data /tmp/test.ffindex nonexistquery b foo || echo
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" ./ffindex_get  /nonexist/data /tmp/test.ffindex a foo || echo
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" ./ffindex_get  /tmp/test.data /nonexist/index a foo || echo

time:
	LD_LIBRARY_PATH=".:$(LD_LIBRARY_PATH)" zsh -c "time ./ffindex_get  ../test_big/nr20_11Jan10_a3m_db ../test_big/nr20_11Jan10_a3m_db.ffindex BABBAMIBA.a3m KARDUZEBA.a3m HARPAQIBA.a3m WOQQATEBA.a3m BEBRAQIBA.a3m GIRSAKEBA.a3m NAQTIFIBA.a3m BIWCIQABA.a3m > /dev/null"

clean:
	rm -f ffindex.o ffutil.o ffindex_build ffindex_build.o ffindex_get ffindex_get.o libffindex.so libffindex.so.0.1 ffindex_modify ffindex_modify.o \
		libffindex.a  ffindex_from_fasta.o ffindex_from_fasta ffindex_apply_mpi.o ffindex_apply.o

install:
	mkdir -p $(INSTALL_DIR)/bin
	mkdir -p $(INSTALL_DIR)/$(libdir)
	mkdir -p $(INSTALL_DIR)/include
	install ffindex_get $(INSTALL_DIR)/bin/ffindex_get
	install ffindex_build $(INSTALL_DIR)/bin/ffindex_build
	install ffindex_modify $(INSTALL_DIR)/bin/ffindex_modify
	install ffindex_from_fasta $(INSTALL_DIR)/bin/ffindex_from_fasta
	install ffindex_apply $(INSTALL_DIR)/bin/ffindex_apply
	install ffindex_unpack $(INSTALL_DIR)/bin/ffindex_unpack
	install libffindex.so.0.1 $(INSTALL_DIR)/$(libdir)/libffindex.so.0.1
	install libffindex.a $(INSTALL_DIR)/$(libdir)/libffindex.a
	ln -sf libffindex.so.0.1 $(INSTALL_DIR)/$(libdir)/libffindex.so
	install ffindex.h $(INSTALL_DIR)/include/ffindex.h
	install ffutil.h $(INSTALL_DIR)/include/ffutil.h
ifdef HAVE_MPI
	install ffindex_apply_mpi $(INSTALL_DIR)/bin/ffindex_apply_mpi
endif

deinstall:
	rm -f $(INSTALL_DIR)/bin/ffindex_get $(INSTALL_DIR)/bin/ffindex_build $(INSTALL_DIR)/bin/ffindex_modify \
		$(INSTALL_DIR)/$(libdir)/libffindex.so.0.1 $(INSTALL_DIR)/$(libdir)/libffindex.a \
		$(INSTALL_DIR)/$(libdir)/libffindex.so $(INSTALL_DIR)/include/ffindex.h $(INSTALL_DIR)/include/ffutil.h\
		$(INSTALL_DIR)/bin/ffindex_from_fasta
	rmdir $(INSTALL_DIR)/bin || true
	rmdir $(INSTALL_DIR)/$(libdir) || true
	rmdir $(INSTALL_DIR)/include || true
