OS:= $(shell uname)

ifeq ($(OS), Darwin)
MFILE=Makefile.osx
else
MFILE=Makefile
endif

all:
	$(MAKE) -C src -f $(MFILE) $@

%:
	$(MAKE) -C src -f $(MFILE) $@

release:
	hg archive -t tgz ffindex-`cat VERSION`.tar.gz
