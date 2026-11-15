#ifndef LIBMARV_AMD_ARCH_CUH
#define LIBMARV_AMD_ARCH_CUH

#if !defined(__CUDACC__)

#include <cstring>
#include <string>

#include "../cuda_errorcheck.cuh"
#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

// Tuned AMD config tables (aligners/make_aligner_amd_<table>.cu)
enum class AmdConfigTable{ none, gfx900, gfx90a, gfx942, gfx950, gfx1100, gfx1150 };

// Table for the device by gcnArchName; must match marv_hip_tables_for_arch in cmake/MarvGpuArch.cmake
inline AmdConfigTable amdConfigTable(int deviceId){
    cudaDeviceProp prop;
    CUDACHECK(cudaGetDeviceProperties(&prop, deviceId));
    const std::string arch(prop.gcnArchName, std::strcspn(prop.gcnArchName, ":"));
    auto startsWith = [&](const char* prefix){ return arch.rfind(prefix, 0) == 0; };

    if(startsWith("gfx115")) return AmdConfigTable::gfx1150;
    if(startsWith("gfx10") || startsWith("gfx11") || startsWith("gfx12")) return AmdConfigTable::gfx1100;
    if(arch == "gfx90a") return AmdConfigTable::gfx90a;
    if(startsWith("gfx95")) return AmdConfigTable::gfx950;
    if(startsWith("gfx94")) return AmdConfigTable::gfx942;
    if(startsWith("gfx9")) return AmdConfigTable::gfx900;
    return AmdConfigTable::none;
}

LIBMARV_NAMESPACE_END

#endif

#endif
