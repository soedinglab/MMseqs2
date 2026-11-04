
NVCC ?= nvcc

#to compile for multiple archs, specify space separated list of archs
#e.g. make GPUARCH="89 90 120"

GPUARCH ?= 120

ARCHS_THAT_HAVE_FAMILY = 100 103 110 120 121

MAKE_ARCH_ARG_NO_FAMILY = -gencode=arch=compute_$1,code=sm_$1
MAKE_ARCH_ARG_WITH_FAMILY = -gencode=arch=compute_$1f,code=sm_$1f



NVCC_GPU_ARCH := $(foreach N, $(GPUARCH), $(if $(findstring $(N), $(ARCHS_THAT_HAVE_FAMILY)), $(call MAKE_ARCH_ARG_WITH_FAMILY,$(N)), $(call MAKE_ARCH_ARG_NO_FAMILY,$(N))))



GPUARCH_NUM_COMPILE_THREADS ?= 1
COMMON_FLAGS = -O3 -std=c++17 --expt-relaxed-constexpr --extended-lambda -Xcompiler="-fopenmp" #-lineinfo -g
COMMON_FLAGS += --threads $(GPUARCH_NUM_COMPILE_THREADS)

SCORE_ONLY_OBJS := benchmark_uint8x4_scoreonly_singleTile.o \
	benchmark_uint8x4_scoreonly_multiTile.o \
	benchmark_short2_scoreonly_singleTile.o \
	benchmark_short2_scoreonly_multiTile.o \
	benchmark_half2_scoreonly_singleTile.o \
	benchmark_half2_scoreonly_multiTile.o \
	benchmark_float_scoreonly_singleTile.o \
	benchmark_float_scoreonly_multiTile.o \
	benchmark_int_scoreonly_singleTile.o \
	benchmark_int_scoreonly_multiTile.o

SCORE_ENDPOS_OBJS := benchmark_uint8x4_scoreonly_singleTile.o \
	benchmark_uint8x4_scoreendpos_multiTile.o \
	benchmark_short2_scoreendpos_singleTile.o \
	benchmark_short2_scoreendpos_multiTile.o \
	benchmark_half2_scoreendpos_singleTile.o \
	benchmark_half2_scoreendpos_multiTile.o \
	benchmark_float_scoreendpos_singleTile.o \
	benchmark_float_scoreendpos_multiTile.o \
	benchmark_int_scoreendpos_singleTile.o \
	benchmark_int_scoreendpos_multiTile.o


tuning_gapless: tuning_gapless.cu benchmarkdata.o $(SCORE_ONLY_OBJS) $(SCORE_ENDPOS_OBJS)
	$(NVCC) $(NVCC_GPU_ARCH) $(COMMON_FLAGS) $^ -o $@


benchmarkdata.o: benchmarkdata.cu benchmarkdata.hpp common.hpp
	nvcc benchmarkdata.cu $(NVCC_GPU_ARCH) $(COMMON_FLAGS) -c -o $@

$(SCORE_ONLY_OBJS): %.o: %.cu benchmark_scoreonly.cuh benchmarkdata.hpp common.hpp
	$(NVCC) $(NVCC_GPU_ARCH) $(COMMON_FLAGS) -c $< -o $@

$(SCORE_ENDPOS_OBJS): %.o: %.cu benchmark_scoreendpos.cuh benchmarkdata.hpp common.hpp
	$(NVCC) $(NVCC_GPU_ARCH) $(COMMON_FLAGS) -c $< -o $@


.PHONY: clean 
clean:
	rm -rf *.o tuning_gapless
	