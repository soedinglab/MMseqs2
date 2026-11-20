# settings
DIALECT      = -std=c++17
OPTIMIZATION = -O3 -g
WARNINGS     = -Xcompiler="-Wall -Wextra"

# NVCC_FLAGS   = -DCUDASW_DEBUG_CHECK_CORRECTNESS -arch=native -lineinfo --expt-relaxed-constexpr -rdc=true --extended-lambda -lnvToolsExt -Xcompiler="-fopenmp" #-res-usage #-Xptxas "-v"

GPUARCH ?= native
NVCC_ARCH = -arch=native

MARVDEFS = -DMARV_HAVE_NVIDIA_SM120

ifeq ($(GPUARCH), native)
	NVCC_ARCH = -arch=native
else
	NVCC_ARCH := $(foreach N, $(GPUARCH),-gencode=arch=compute_$(N),code=sm_$(N) )
	MARVDEFS := $(foreach N, $(GPUARCH),-DMARV_HAVE_NVIDIA_SM$(N))
endif

$(info MARVDEFS: $(MARVDEFS))

NVCC_FLAGS   =  $(MARVDEFS) -lineinfo --expt-relaxed-constexpr --extended-lambda -Xcompiler="-fopenmp" #-DCUDASW_DEBUG_CHECK_CORRECTNESS # -rdc=true  #-res-usage #-Xptxas "-v"


LDFLAGS      = -Xcompiler="-pthread"  $(NVCC_FLAGS) -lz
NVCC ?= nvcc
COMPILER     = $(NVCC)
ARTIFACT     = align

BUILDDIR = build

MAKEDB = makedb
MODIFYDB = modifydb
GRIDSEARCH = gridsearch
TILECONFIGSEARCH = tileconfigsearch

$(shell mkdir -p $(BUILDDIR))

# make targets
.PHONY: clean

release: $(ARTIFACT) $(MAKEDB)


clean :
	rm -f $(BUILDDIR)/*
	rm -f $(ARTIFACT)
	rm -f $(MAKEDB)

# compiler call
COMPILE = $(COMPILER) $(NVCC_FLAGS) $(NVCC_ARCH) $(DIALECT) $(OPTIMIZATION) $(WARNINGS) -c $< -o $@
COMPILE_MANUAL_ARCH = $(COMPILER) $(NVCC_FLAGS) $(DIALECT) $(OPTIMIZATION) $(WARNINGS) -c $< -o $@

LIBMARV_OBJS =  $(BUILDDIR)/sequence_io.o $(BUILDDIR)/dbdata.o $(BUILDDIR)/options.o $(BUILDDIR)/gpudatabase.o $(BUILDDIR)/dbbatching.o \
	$(BUILDDIR)/score_only_search.o $(BUILDDIR)/score_endpos_search.o $(BUILDDIR)/score_only_foldseek_search.o


LIBMARV_OBJS += $(BUILDDIR)/make_aligner.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm75.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm80.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm86.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm87.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm89.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm90.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm100.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm103.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm110.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm120.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_sm121.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_defaultptx.o

LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm75.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm80.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm86.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm87.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm89.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm90.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm100.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm103.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm110.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm120.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_sm121.o
LIBMARV_OBJS += $(BUILDDIR)/make_aligner_endpos_defaultptx.o


# link object files into executable
$(ARTIFACT): $(BUILDDIR)/main2.o $(LIBMARV_OBJS)
	$(COMPILER) $^ -o $(ARTIFACT) $(LDFLAGS)

$(MAKEDB): $(BUILDDIR)/makedb.o $(BUILDDIR)/sequence_io.o $(BUILDDIR)/dbdata.o
	$(COMPILER) $^ -o $(MAKEDB) $(LDFLAGS)


main2: $(BUILDDIR)/main2.o $(LIBMARV_OBJS)
	$(COMPILER) $^ -o main2 $(LDFLAGS)

$(BUILDDIR)/main2.o : src/main2.cu src/gpudatabase.cuh src/searches/score_only_search.cuh src/searches/score_endpos_search.cuh src/searches/score_only_foldseek_search.cuh
	$(COMPILE)



$(BUILDDIR)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h
	$(COMPILE)

$(BUILDDIR)/dbdata.o : src/dbdata.cpp src/dbdata.hpp src/mapped_file.hpp src/sequence_io.h src/length_partitions.hpp
	$(COMPILE)

$(BUILDDIR)/options.o : src/options.cpp src/options.hpp src/types.hpp
	$(COMPILE)

$(BUILDDIR)/makedb.o : src/makedb.cpp src/dbdata.hpp src/sequence_io.h
	$(COMPILE)

$(BUILDDIR)/marv.o : src/marv.cu src/marv.h src/dbdata.hpp src/gpudatabase.cuh src/config.hpp src/base64.h
	$(COMPILE)

$(BUILDDIR)/gpudatabase.o : src/gpudatabase.cu src/gpudatabase.cuh
	$(COMPILE)

$(BUILDDIR)/dbbatching.o : src/dbbatching.cu src/dbbatching.cuh
	$(COMPILE)

$(BUILDDIR)/score_only_search.o : src/searches/score_only_search.cu src/searches/score_only_search.cuh
	$(COMPILE)

$(BUILDDIR)/score_endpos_search.o : src/searches/score_endpos_search.cu src/searches/score_endpos_search.cuh
	$(COMPILE)

$(BUILDDIR)/score_only_foldseek_search.o : src/searches/score_only_foldseek_search.cu src/searches/score_only_foldseek_search.cuh
	$(COMPILE)


$(BUILDDIR)/make_aligner_sm75.o : src/aligners/make_aligner_sm75.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_75,code=sm_75

$(BUILDDIR)/make_aligner_sm80.o : src/aligners/make_aligner_sm80.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_80,code=sm_80

$(BUILDDIR)/make_aligner_sm86.o : src/aligners/make_aligner_sm86.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_86,code=sm_86

$(BUILDDIR)/make_aligner_sm87.o : src/aligners/make_aligner_sm87.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_87,code=sm_87

$(BUILDDIR)/make_aligner_sm89.o : src/aligners/make_aligner_sm89.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_89,code=sm_89

$(BUILDDIR)/make_aligner_sm90.o : src/aligners/make_aligner_sm90.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_90,code=sm_90

$(BUILDDIR)/make_aligner_sm100.o : src/aligners/make_aligner_sm100.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_100,code=sm_100

$(BUILDDIR)/make_aligner_sm103.o : src/aligners/make_aligner_sm103.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_103,code=sm_103

$(BUILDDIR)/make_aligner_sm110.o : src/aligners/make_aligner_sm110.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_110,code=sm_110

$(BUILDDIR)/make_aligner_sm120.o : src/aligners/make_aligner_sm120.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_120f,code=sm_120

$(BUILDDIR)/make_aligner_sm121.o : src/aligners/make_aligner_sm121.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_121f,code=sm_121

$(BUILDDIR)/make_aligner_defaultptx.o : src/aligners/make_aligner_defaultptx.cu src/aligners/score_only_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_80,code=compute_80

$(BUILDDIR)/make_aligner.o : src/aligners/make_aligner.cu
	$(COMPILE)



$(BUILDDIR)/make_aligner_endpos_sm75.o : src/aligners/make_aligner_endpos_sm75.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_75,code=sm_75

$(BUILDDIR)/make_aligner_endpos_sm80.o : src/aligners/make_aligner_endpos_sm80.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_80,code=sm_80

$(BUILDDIR)/make_aligner_endpos_sm86.o : src/aligners/make_aligner_endpos_sm86.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_86,code=sm_86

$(BUILDDIR)/make_aligner_endpos_sm87.o : src/aligners/make_aligner_endpos_sm87.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_87,code=sm_87

$(BUILDDIR)/make_aligner_endpos_sm89.o : src/aligners/make_aligner_endpos_sm89.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_89,code=sm_89

$(BUILDDIR)/make_aligner_endpos_sm90.o : src/aligners/make_aligner_endpos_sm90.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_90,code=sm_90

$(BUILDDIR)/make_aligner_endpos_sm100.o : src/aligners/make_aligner_endpos_sm100.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_100,code=sm_100

$(BUILDDIR)/make_aligner_endpos_sm103.o : src/aligners/make_aligner_endpos_sm103.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_103,code=sm_103

$(BUILDDIR)/make_aligner_endpos_sm110.o : src/aligners/make_aligner_endpos_sm110.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_110,code=sm_110

$(BUILDDIR)/make_aligner_endpos_sm120.o : src/aligners/make_aligner_endpos_sm120.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_120f,code=sm_120

$(BUILDDIR)/make_aligner_endpos_sm121.o : src/aligners/make_aligner_endpos_sm121.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_121f,code=sm_121

$(BUILDDIR)/make_aligner_endpos_defaultptx.o : src/aligners/make_aligner_endpos_defaultptx.cu src/aligners/score_with_endpos_alignment.cuh
	$(COMPILE_MANUAL_ARCH) -gencode=arch=compute_80,code=compute_80

$(BUILDDIR)/make_aligner_endpos.o : src/aligners/make_aligner_endpos.cu
	$(COMPILE)


	



