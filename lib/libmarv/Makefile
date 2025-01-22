# settings
DIALECT      = -std=c++17
OPTIMIZATION = -O3 -g
WARNINGS     = -Xcompiler="-Wall -Wextra"

# NVCC_FLAGS   = -DCUDASW_DEBUG_CHECK_CORRECTNESS -arch=native -lineinfo --expt-relaxed-constexpr -rdc=true --extended-lambda -lnvToolsExt -Xcompiler="-fopenmp" #-res-usage #-Xptxas "-v"
NVCC_FLAGS   = -arch=native -lineinfo --expt-relaxed-constexpr -rdc=true --extended-lambda -lnvToolsExt -Xcompiler="-fopenmp" #-res-usage #-Xptxas "-v"


LDFLAGS      = -Xcompiler="-pthread"  $(NVCC_FLAGS) -lz
COMPILER     = nvcc
ARTIFACT     = align

BUILDDIR = build

MAKEDB = makedb
MODIFYDB = modifydb
GRIDSEARCH = gridsearch
TILECONFIGSEARCH = tileconfigsearch

$(shell mkdir -p $(BUILDDIR))

# make targets
.PHONY: clean

release: $(ARTIFACT) $(MAKEDB) $(TILECONFIGSEARCH)


clean :
	rm -f $(BUILDDIR)/*
	rm -f $(ARTIFACT)
	rm -f $(MAKEDB)
	rm -f $(TILECONFIGSEARCH)

# compiler call
COMPILE = $(COMPILER) $(NVCC_FLAGS) $(DIALECT) $(OPTIMIZATION) $(WARNINGS) -c $< -o $@

CUDASW_OBJS =  $(BUILDDIR)/sequence_io.o $(BUILDDIR)/dbdata.o $(BUILDDIR)/options.o $(BUILDDIR)/blosum.o $(BUILDDIR)/pssmkernels_smithwaterman_instantiation_float.o $(BUILDDIR)/pssmkernels_smithwaterman_instantiation_dpx.o $(BUILDDIR)/pssmkernels_gapless_instantiation_half2.o $(BUILDDIR)/pssmkernels_gapless_instantiation_dpx.o $(BUILDDIR)/pssmkernels_gapless_instantiation_half2_kernelparamzero.o $(BUILDDIR)/pssmkernels_gapless_instantiation_dpx_kernelparamzero.o

# link object files into executable
$(ARTIFACT): $(BUILDDIR)/main.o $(CUDASW_OBJS)
	$(COMPILER) $^ -o $(ARTIFACT) $(LDFLAGS)

$(TILECONFIGSEARCH): $(BUILDDIR)/tileconfigsearch.o $(CUDASW_OBJS)
	$(COMPILER) $^ -o $(TILECONFIGSEARCH) $(LDFLAGS)	

$(MAKEDB): $(BUILDDIR)/makedb.o $(BUILDDIR)/sequence_io.o $(BUILDDIR)/dbdata.o
	$(COMPILER) $^ -o $(MAKEDB) $(LDFLAGS)



$(BUILDDIR)/tileconfigsearch.o : src/tileconfigsearch.cu src/sequence_io.h src/length_partitions.hpp src/dbdata.hpp src/cudasw4.cuh src/kernels.cuh src/convert.cuh src/blosum.hpp src/types.hpp src/pssm.cuh src/pssmkernels_gapless.cuh src/pssmkernels_smithwaterman.cuh src/gapless_kernel_config.cuh
	$(COMPILE)

$(BUILDDIR)/main.o : src/main.cu src/sequence_io.h src/length_partitions.hpp src/dbdata.hpp src/cudasw4.cuh src/kernels.cuh src/convert.cuh src/blosum.hpp src/types.hpp src/pssm.cuh src/pssmkernels_gapless.cuh src/pssmkernels_smithwaterman.cuh src/gapless_kernel_config.cuh
	$(COMPILE)

$(BUILDDIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h
	$(COMPILE)

$(BUILDDIR)/dbdata.o : src/dbdata.cpp src/dbdata.hpp src/mapped_file.hpp src/sequence_io.h src/length_partitions.hpp
	$(COMPILE)

$(BUILDDIR)/options.o : src/options.cpp src/options.hpp src/types.hpp
	$(COMPILE)

$(BUILDDIR)/blosum.o : src/blosum.cu src/blosum.hpp
	$(COMPILE)

$(BUILDDIR)/pssmkernels_gapless_instantiation_half2.o : src/pssmkernels_gapless_instantiation_half2.cu src/pssmkernels_gapless.cuh src/convert.cuh src/util.cuh
	$(COMPILE)

$(BUILDDIR)/pssmkernels_gapless_instantiation_dpx.o : src/pssmkernels_gapless_instantiation_dpx.cu src/pssmkernels_gapless.cuh src/convert.cuh src/util.cuh
	$(COMPILE)	

$(BUILDDIR)/pssmkernels_gapless_instantiation_half2_kernelparamzero.o : src/pssmkernels_gapless_instantiation_half2_kernelparamzero.cu src/pssmkernels_gapless.cuh src/convert.cuh src/util.cuh
	$(COMPILE)

$(BUILDDIR)/pssmkernels_gapless_instantiation_dpx_kernelparamzero.o : src/pssmkernels_gapless_instantiation_dpx_kernelparamzero.cu src/pssmkernels_gapless.cuh src/convert.cuh src/util.cuh
	$(COMPILE)	

$(BUILDDIR)/pssmkernels_smithwaterman_instantiation_float.o : src/pssmkernels_smithwaterman_instantiation_float.cu src/pssmkernels_smithwaterman.cuh src/convert.cuh src/util.cuh
	$(COMPILE)

$(BUILDDIR)/pssmkernels_smithwaterman_instantiation_dpx.o : src/pssmkernels_smithwaterman_instantiation_dpx.cu src/pssmkernels_smithwaterman.cuh src/convert.cuh src/util.cuh
	$(COMPILE)

$(BUILDDIR)/makedb.o : src/makedb.cpp src/dbdata.hpp src/sequence_io.h
	$(COMPILE)



