#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>


#include "make_aligner_headers.cuh"
#include "../cuda_errorcheck.cuh"

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_32bit(){
    int deviceId = 0;
    CUDACHECK(cudaGetDevice(&deviceId));
    int ccMajor = 0;
    int ccMinor = 0;
    CUDACHECK(cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&ccMinor, cudaDevAttrComputeCapabilityMinor, deviceId));

#if defined(__CUDACC__)

#ifdef MARV_HAVE_NVIDIA_SM75
    if(ccMajor == 7 && ccMinor == 5){
        return make_multiconfig_sw_endpos_aligner_sm75_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM80
    if(ccMajor == 8 && ccMinor == 0){
        return make_multiconfig_sw_endpos_aligner_sm80_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM86
    if(ccMajor == 8 && ccMinor == 6){
        return make_multiconfig_sw_endpos_aligner_sm86_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM87
    if(ccMajor == 8 && ccMinor == 7){
        return make_multiconfig_sw_endpos_aligner_sm87_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM89
    if(ccMajor == 8 && ccMinor == 9){
        return make_multiconfig_sw_endpos_aligner_sm89_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM90
    if(ccMajor == 9 && ccMinor == 0){
        return make_multiconfig_sw_endpos_aligner_sm90_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM100
    if(ccMajor == 10 && ccMinor == 0){
        return make_multiconfig_sw_endpos_aligner_sm100_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM103
    if(ccMajor == 10 && ccMinor == 3){
        return make_multiconfig_sw_endpos_aligner_sm103_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM110
    if((ccMajor == 11 && ccMinor == 0) || (ccMajor == 10 && ccMinor == 1)){
        return make_multiconfig_sw_endpos_aligner_sm110_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM120
    if(ccMajor == 12 && ccMinor == 0){
        return make_multiconfig_sw_endpos_aligner_sm120_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM121
    if(ccMajor == 12 && ccMinor == 1){
        return make_multiconfig_sw_endpos_aligner_sm121_32bit();
    }
#endif
    return make_multiconfig_sw_endpos_aligner_default_ptx_32bit();

#else


#ifdef MARV_HAVE_AMD_GFX900
    if(ccMajor == 9 && ccMinor == 0){
        return make_multiconfig_sw_endpos_aligner_amd_gfx900_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX942 
    if(ccMajor == 9 && ccMinor == 4){
        return make_multiconfig_sw_endpos_aligner_amd_gfx942_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX1100 
    if(ccMajor == 11 && ccMinor == 0){
        return make_multiconfig_sw_endpos_aligner_amd_gfx1100_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX1150
    if(ccMajor == 11 && ccMinor == 5){
        return make_multiconfig_sw_endpos_aligner_amd_gfx1150_32bit();
    }
#endif

    {
        //return make_multiconfig_sw_endpos_aligner_default_amd_32bit();
        throw std::runtime_error("No aligner found for detected AMD arch: major " 
                        + std::to_string(ccMajor) + ", minor " + std::to_string(ccMinor));
    }

    

#endif
}



std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_16bit(){
    int deviceId = 0;
    CUDACHECK(cudaGetDevice(&deviceId));
    int ccMajor = 0;
    int ccMinor = 0;
    CUDACHECK(cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&ccMinor, cudaDevAttrComputeCapabilityMinor, deviceId));

#if defined(__CUDACC__)

#ifdef MARV_HAVE_NVIDIA_SM75
    if(ccMajor == 7 && ccMinor == 5){
        return make_multiconfig_gapless_endpos_aligner_sm75_32bit(); //does not support 16bit
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM80
    if(ccMajor == 8 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_sm80_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM86
    if(ccMajor == 8 && ccMinor == 6){
        return make_multiconfig_gapless_endpos_aligner_sm86_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM87
    if(ccMajor == 8 && ccMinor == 7){
        return make_multiconfig_gapless_endpos_aligner_sm87_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM89
    if(ccMajor == 8 && ccMinor == 9){
        return make_multiconfig_gapless_endpos_aligner_sm89_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM90
    if(ccMajor == 9 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_sm90_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM100
    if(ccMajor == 10 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_sm100_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM103
    if(ccMajor == 10 && ccMinor == 3){
        return make_multiconfig_gapless_endpos_aligner_sm103_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM110
    if((ccMajor == 11 && ccMinor == 0) || (ccMajor == 10 && ccMinor == 1)){
        return make_multiconfig_gapless_endpos_aligner_sm110_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM120
    if(ccMajor == 12 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_sm120_16bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM121
    if(ccMajor == 12 && ccMinor == 1){
        return make_multiconfig_gapless_endpos_aligner_sm121_16bit();
    }
#endif

    return make_multiconfig_gapless_endpos_aligner_default_ptx_16bit();

#else


#ifdef MARV_HAVE_AMD_GFX900
    if(ccMajor == 9 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_amd_gfx900_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX942
    if(ccMajor == 9 && ccMinor == 4){
        return make_multiconfig_gapless_endpos_aligner_amd_gfx942_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX1100
    if(ccMajor == 11 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_amd_gfx1100_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX1150
    if(ccMajor == 11 && ccMinor == 5){
        return make_multiconfig_gapless_endpos_aligner_amd_gfx1150_32bit();
    }
#endif

    {
        // return make_multiconfig_gapless_endpos_aligner_default_amd_32bit();
        throw std::runtime_error("No aligner found for detected AMD arch: major " 
                        + std::to_string(ccMajor) + ", minor " + std::to_string(ccMinor));
    }

#endif
}

std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_32bit(){
    int deviceId = 0;
    CUDACHECK(cudaGetDevice(&deviceId));
    int ccMajor = 0;
    int ccMinor = 0;
    CUDACHECK(cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&ccMinor, cudaDevAttrComputeCapabilityMinor, deviceId));

#if defined(__CUDACC__)

#ifdef MARV_HAVE_NVIDIA_SM75
    if(ccMajor == 7 && ccMinor == 5){
        return make_multiconfig_gapless_endpos_aligner_sm75_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM80
    if(ccMajor == 8 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_sm80_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM86
    if(ccMajor == 8 && ccMinor == 6){
        return make_multiconfig_gapless_endpos_aligner_sm86_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM87
    if(ccMajor == 8 && ccMinor == 7){
        return make_multiconfig_gapless_endpos_aligner_sm87_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM89
    if(ccMajor == 8 && ccMinor == 9){
        return make_multiconfig_gapless_endpos_aligner_sm89_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM90
    if(ccMajor == 9 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_sm90_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM100
    if(ccMajor == 10 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_sm100_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM103
    if(ccMajor == 10 && ccMinor == 3){
        return make_multiconfig_gapless_endpos_aligner_sm103_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM110
    if((ccMajor == 11 && ccMinor == 0) || (ccMajor == 10 && ccMinor == 1)){
        return make_multiconfig_gapless_endpos_aligner_sm110_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM120
    if(ccMajor == 12 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_sm120_32bit();
    }
#endif
#ifdef MARV_HAVE_NVIDIA_SM121
    if(ccMajor == 12 && ccMinor == 1){
        return make_multiconfig_gapless_endpos_aligner_sm121_32bit();
    }
#endif

    return make_multiconfig_gapless_endpos_aligner_default_ptx_32bit();

#else


#ifdef MARV_HAVE_AMD_GFX900
    if(ccMajor == 9 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_amd_gfx900_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX942
    if(ccMajor == 9 && ccMinor == 4){
        return make_multiconfig_gapless_endpos_aligner_amd_gfx942_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX1100
    if(ccMajor == 11 && ccMinor == 0){
        return make_multiconfig_gapless_endpos_aligner_amd_gfx1100_32bit();
    }
#endif
#ifdef MARV_HAVE_AMD_GFX1150
    if(ccMajor == 11 && ccMinor == 5){
        return make_multiconfig_gapless_endpos_aligner_amd_gfx1150_32bit();
    }
#endif

    {
        // return make_multiconfig_gapless_endpos_aligner_default_amd_32bit();
        throw std::runtime_error("No aligner found for detected AMD arch: major " 
                        + std::to_string(ccMajor) + ", minor " + std::to_string(ccMinor));
    }

#endif


}

LIBMARV_NAMESPACE_END