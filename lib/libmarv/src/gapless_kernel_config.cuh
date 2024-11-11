#ifndef GAPLESS_KERNEL_CONFIG_CUH
#define GAPLESS_KERNEL_CONFIG_CUH

#include <algorithm>
#include <vector>
#include <string>

namespace cudasw4{


    struct GaplessKernelConfig{
        enum class Approach : int{
            hardcodedzero = 0,
            kernelparamzero = 1
        };
        bool dpx;
        int tilesize;
        int groupsize;
        int numRegs;
        Approach approach;

        GaplessKernelConfig() = default;
        GaplessKernelConfig(int tilesize_, int groupsize_, int numRegs_, int dpx_, Approach approach_)
            : dpx(dpx_), tilesize(tilesize_), groupsize(groupsize_), numRegs(numRegs_), approach(approach_)
        {}

        GaplessKernelConfig(const GaplessKernelConfig&) = default;
        GaplessKernelConfig& operator=(const GaplessKernelConfig&) = default;
    };

    __inline__
    std::string to_string(GaplessKernelConfig::Approach approach){
        switch(approach){
            case GaplessKernelConfig::Approach::hardcodedzero: return "hardcodedzero";
            case GaplessKernelConfig::Approach::kernelparamzero: return "kernelparamzero";
        }
        return "to_string: missing case for GaplessKernelConfig::Approach";
    }

    __inline__
    std::ostream& operator<<(std::ostream& os, const GaplessKernelConfig& data){

        os << data.tilesize << " " << data.groupsize << " " << data.numRegs 
            << " " << data.dpx << " " << int(data.approach);
        return os;
    }
    

    //T4
    __inline__
    std::vector<GaplessKernelConfig> getOptimalKernelConfigs_gapless_sm75(){
        std::vector<GaplessKernelConfig> configs{
            {32,4,4,0,GaplessKernelConfig::Approach::hardcodedzero},
            {64,4,8,0,GaplessKernelConfig::Approach::hardcodedzero},
            {96,4,12,0,GaplessKernelConfig::Approach::hardcodedzero},
            {128,4,16,0,GaplessKernelConfig::Approach::hardcodedzero},
            {160,4,20,0,GaplessKernelConfig::Approach::hardcodedzero},
            {192,8,12,0,GaplessKernelConfig::Approach::hardcodedzero},
            {224,4,28,0,GaplessKernelConfig::Approach::hardcodedzero},
            {256,8,16,0,GaplessKernelConfig::Approach::hardcodedzero},
            {288,4,36,0,GaplessKernelConfig::Approach::hardcodedzero},
            {320,8,20,0,GaplessKernelConfig::Approach::hardcodedzero},
            {352,4,44,0,GaplessKernelConfig::Approach::hardcodedzero},
            {384,16,12,0,GaplessKernelConfig::Approach::hardcodedzero},
            {416,4,52,0,GaplessKernelConfig::Approach::hardcodedzero},
            {448,8,28,0,GaplessKernelConfig::Approach::hardcodedzero},
            {480,4,60,0,GaplessKernelConfig::Approach::hardcodedzero},
            {512,16,16,0,GaplessKernelConfig::Approach::hardcodedzero},
            {576,8,36,0,GaplessKernelConfig::Approach::hardcodedzero},
            {640,16,20,0,GaplessKernelConfig::Approach::hardcodedzero},
            {704,8,44,0,GaplessKernelConfig::Approach::hardcodedzero},
            {768,16,24,0,GaplessKernelConfig::Approach::hardcodedzero},
            {832,8,52,0,GaplessKernelConfig::Approach::hardcodedzero},
            {896,16,28,0,GaplessKernelConfig::Approach::hardcodedzero},
            {960,8,60,0,GaplessKernelConfig::Approach::hardcodedzero},
            {1024,16,32,0,GaplessKernelConfig::Approach::hardcodedzero},
            {1152,16,36,0,GaplessKernelConfig::Approach::hardcodedzero},
            {1280,16,40,0,GaplessKernelConfig::Approach::hardcodedzero},
            {1408,16,44,0,GaplessKernelConfig::Approach::hardcodedzero},
            {1536,16,48,0,GaplessKernelConfig::Approach::hardcodedzero},
            //larger tiles are not supported because shared memory size is too small
        };

        return configs;
    }


    //A100
    __inline__
    std::vector<GaplessKernelConfig> getOptimalKernelConfigs_gapless_sm80(){
        std::vector<GaplessKernelConfig> configs{
            {32,4,4,0, GaplessKernelConfig::Approach::kernelparamzero},
            {64,4,8,0, GaplessKernelConfig::Approach::hardcodedzero},
            {96,4,12,0, GaplessKernelConfig::Approach::kernelparamzero},
            {128,4,16,0, GaplessKernelConfig::Approach::kernelparamzero},
            {160,4,20,0, GaplessKernelConfig::Approach::kernelparamzero},
            {192,4,24,0, GaplessKernelConfig::Approach::kernelparamzero},
            {224,4,28,0, GaplessKernelConfig::Approach::kernelparamzero},
            {256,4,32,0, GaplessKernelConfig::Approach::kernelparamzero},
            {288,4,36,0, GaplessKernelConfig::Approach::kernelparamzero},
            {320,4,40,0, GaplessKernelConfig::Approach::kernelparamzero},
            {352,4,44,0, GaplessKernelConfig::Approach::kernelparamzero},
            {384,4,48,0, GaplessKernelConfig::Approach::kernelparamzero},
            {416,4,52,0, GaplessKernelConfig::Approach::kernelparamzero},
            {448,4,56,0, GaplessKernelConfig::Approach::hardcodedzero},
            {480,4,60,0, GaplessKernelConfig::Approach::hardcodedzero},
            {512,4,64,0, GaplessKernelConfig::Approach::hardcodedzero},
            {576,8,36,0, GaplessKernelConfig::Approach::kernelparamzero},
            {640,8,40,0, GaplessKernelConfig::Approach::hardcodedzero},
            {704,8,44,0, GaplessKernelConfig::Approach::hardcodedzero},
            {768,8,48,0, GaplessKernelConfig::Approach::kernelparamzero},
            {832,8,52,0, GaplessKernelConfig::Approach::hardcodedzero},
            {896,8,56,0, GaplessKernelConfig::Approach::kernelparamzero},
            {960,8,60,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1024,8,64,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1152,16,36,0, GaplessKernelConfig::Approach::kernelparamzero},
            {1280,16,40,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1408,16,44,0, GaplessKernelConfig::Approach::kernelparamzero},
            {1536,16,48,0, GaplessKernelConfig::Approach::kernelparamzero},
            {1664,16,52,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1792,16,56,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1920,16,60,0, GaplessKernelConfig::Approach::kernelparamzero},
            {2048,16,64,0, GaplessKernelConfig::Approach::kernelparamzero},
        };

        return configs;
    }

    //L40S
    __inline__
    std::vector<GaplessKernelConfig> getOptimalKernelConfigs_gapless_sm89(){
        std::vector<GaplessKernelConfig> configs{
            {32,4,4,0, GaplessKernelConfig::Approach::kernelparamzero},
            {64,4,8,0, GaplessKernelConfig::Approach::kernelparamzero},
            {96,4,12,0, GaplessKernelConfig::Approach::kernelparamzero},
            {128,4,16,0, GaplessKernelConfig::Approach::kernelparamzero},
            {160,4,20,0, GaplessKernelConfig::Approach::hardcodedzero},
            {192,4,24,0, GaplessKernelConfig::Approach::kernelparamzero},
            {224,4,28,0, GaplessKernelConfig::Approach::hardcodedzero},
            {256,4,32,0, GaplessKernelConfig::Approach::hardcodedzero},
            {288,4,36,0, GaplessKernelConfig::Approach::hardcodedzero},
            {320,4,40,0, GaplessKernelConfig::Approach::hardcodedzero},
            {352,4,44,0, GaplessKernelConfig::Approach::hardcodedzero},
            {384,4,48,0, GaplessKernelConfig::Approach::kernelparamzero},
            {416,4,52,0, GaplessKernelConfig::Approach::hardcodedzero},
            {448,4,56,0, GaplessKernelConfig::Approach::kernelparamzero},
            {480,4,60,0, GaplessKernelConfig::Approach::kernelparamzero},
            {512,4,64,0, GaplessKernelConfig::Approach::hardcodedzero},
            {576,8,36,0, GaplessKernelConfig::Approach::kernelparamzero},
            {640,8,40,0, GaplessKernelConfig::Approach::kernelparamzero},
            {704,8,44,0, GaplessKernelConfig::Approach::kernelparamzero},
            {768,8,48,0, GaplessKernelConfig::Approach::kernelparamzero},
            {832,8,52,0, GaplessKernelConfig::Approach::kernelparamzero},
            {896,8,56,0, GaplessKernelConfig::Approach::kernelparamzero},
            {960,8,60,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1024,8,64,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1152,16,36,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1280,16,40,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1408,16,44,0, GaplessKernelConfig::Approach::kernelparamzero},
            {1536,16,48,0, GaplessKernelConfig::Approach::kernelparamzero},
            {1664,16,52,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1792,16,56,0, GaplessKernelConfig::Approach::hardcodedzero},
            {1920,16,60,0, GaplessKernelConfig::Approach::kernelparamzero},
            {2048,16,64,0, GaplessKernelConfig::Approach::hardcodedzero},
        };

        return configs;
    }

    //H100 SXM
    __inline__
    std::vector<GaplessKernelConfig> getOptimalKernelConfigs_gapless_sm90(){
        std::vector<GaplessKernelConfig> configs{
            {32, 4, 4, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {64, 4, 8, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {96, 4, 12, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {128, 4, 16, 0, GaplessKernelConfig::Approach::kernelparamzero},
            {160, 4, 20, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {192, 4, 24, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {224, 4, 28, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {256, 4, 32, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {288, 4, 36, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {320, 4, 40, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {352, 4, 44, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {384, 4, 48, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {416, 4, 52, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {448, 4, 56, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {480, 4, 60, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {512, 8, 32, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {576, 8, 36, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {640, 8, 40, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {704, 8, 44, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {768, 8, 48, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {832, 8, 52, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {896, 8, 56, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {960, 8, 60, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {1024, 16, 32, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {1152, 16, 36, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {1280, 16, 40, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {1408, 16, 44, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {1536, 16, 48, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {1664, 16, 52, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {1792, 16, 56, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {1920, 16, 60, 1, GaplessKernelConfig::Approach::kernelparamzero},
            {2048, 16, 64, 0, GaplessKernelConfig::Approach::hardcodedzero},
        };

        return configs;
    }

    __inline__
    std::vector<GaplessKernelConfig> getOptimalKernelConfigs_gapless_default(){
        return getOptimalKernelConfigs_gapless_sm89();
    }

    __inline__
    std::vector<GaplessKernelConfig> getOptimalKernelConfigs_gapless(int deviceId){
        int ccMajor = 0;
        int ccMinor = 0;
        cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, deviceId);
        cudaDeviceGetAttribute(&ccMinor, cudaDevAttrComputeCapabilityMinor, deviceId);

        std::vector<GaplessKernelConfig> configs;

        if(ccMajor == 7 && ccMinor == 5){
            configs = getOptimalKernelConfigs_gapless_sm75();
        }else if(ccMajor == 8 && ccMinor == 0){
            configs = getOptimalKernelConfigs_gapless_sm80();
        }else if(ccMajor == 8 && ccMinor == 9){
            configs = getOptimalKernelConfigs_gapless_sm89();
        }else if(ccMajor == 9 && ccMinor == 0){
            configs = getOptimalKernelConfigs_gapless_sm90();
        }else{
            configs = getOptimalKernelConfigs_gapless_default();
        }

        return configs;
    }


}

#endif