#ifndef SW_KERNEL_CONFIG_CUH
#define SW_KERNEL_CONFIG_CUH

#include <algorithm>
#include <vector>
#include <string>

namespace cudasw4{


    struct SmithWatermanKernelConfig{
        enum class Approach : int{
            Unused = 999999
        };
        bool dpx;
        int tilesize;
        int groupsize;
        int numRegs;
        Approach approach;

        SmithWatermanKernelConfig() = default;
        SmithWatermanKernelConfig(int tilesize_, int groupsize_, int numRegs_, int dpx_, Approach approach_)
            : dpx(dpx_), tilesize(tilesize_), groupsize(groupsize_), numRegs(numRegs_), approach(approach_)
        {}

        SmithWatermanKernelConfig(const SmithWatermanKernelConfig&) = default;
        SmithWatermanKernelConfig& operator=(const SmithWatermanKernelConfig&) = default;
    };

    __inline__
    std::string to_string(SmithWatermanKernelConfig::Approach approach){
        switch(approach){
            case SmithWatermanKernelConfig::Approach::Unused: return "Unused";
        }
        return "to_string: missing case for SmithWatermanKernelConfig::Approach";
    }

    __inline__
    std::ostream& operator<<(std::ostream& os, const SmithWatermanKernelConfig& data){

        os << data.tilesize << " " << data.groupsize << " " << data.numRegs 
            << " " << data.dpx << " " << int(data.approach);
        return os;
    }
    
    //T4
    __inline__
    std::vector<SmithWatermanKernelConfig> getOptimalKernelConfigs_SW_sm75(){
        std::vector<SmithWatermanKernelConfig> configs{
            {16,4,4,0,SmithWatermanKernelConfig::Approach::Unused},
            {32,4,8,0,SmithWatermanKernelConfig::Approach::Unused},
            {48,4,12,0,SmithWatermanKernelConfig::Approach::Unused},
            {64,4,16,0,SmithWatermanKernelConfig::Approach::Unused},
            {80,4,20,0,SmithWatermanKernelConfig::Approach::Unused},
            {96,4,24,0,SmithWatermanKernelConfig::Approach::Unused},
            {112,4,28,0,SmithWatermanKernelConfig::Approach::Unused},
            {128,8,16,0,SmithWatermanKernelConfig::Approach::Unused},
            {144,4,36,0,SmithWatermanKernelConfig::Approach::Unused},
            {160,8,20,0,SmithWatermanKernelConfig::Approach::Unused},
            {176,4,44,0,SmithWatermanKernelConfig::Approach::Unused},
            {192,8,24,0,SmithWatermanKernelConfig::Approach::Unused},
            {224,8,28,0,SmithWatermanKernelConfig::Approach::Unused},
            {256,16,16,0,SmithWatermanKernelConfig::Approach::Unused},
            {288,8,36,0,SmithWatermanKernelConfig::Approach::Unused},
            {320,16,20,0,SmithWatermanKernelConfig::Approach::Unused},
            {352,8,44,0,SmithWatermanKernelConfig::Approach::Unused},
            {384,16,24,0,SmithWatermanKernelConfig::Approach::Unused},
            {448,16,28,0,SmithWatermanKernelConfig::Approach::Unused},
            {512,32,16,0,SmithWatermanKernelConfig::Approach::Unused},
            {576,16,36,0,SmithWatermanKernelConfig::Approach::Unused},
            {640,32,20,0,SmithWatermanKernelConfig::Approach::Unused},
            {704,16,44,0,SmithWatermanKernelConfig::Approach::Unused},
            {768,32,24,0,SmithWatermanKernelConfig::Approach::Unused},
            //larger tiles are not supported because shared memory size is too small
        };
        
        return configs;
    }

    //A100
    __inline__
    std::vector<SmithWatermanKernelConfig> getOptimalKernelConfigs_SW_sm80(){
        std::vector<SmithWatermanKernelConfig> configs{
            {16,4,4,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {32,4,8,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {48,4,12,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {64,4,16,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {80,4,20,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {96,4,24,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {112,4,28,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {128,8,16,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {144,4,36,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {160,8,20,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {176,4,44,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {192,8,24,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {224,8,28,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {256,8,32,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {288,8,36,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {320,16,20,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {352,8,44,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {384,16,24,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {448,16,28,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {512,16,32,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {576,16,36,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {640,32,20,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {704,16,44,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {768,32,24,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {896,32,28,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {1024,32,32,0,SmithWatermanKernelConfig::Approach::Unused}, 
        };

        return configs;
    }

    //L40S
    __inline__
    std::vector<SmithWatermanKernelConfig> getOptimalKernelConfigs_SW_sm89(){
        std::vector<SmithWatermanKernelConfig> configs{
            {16,4,4,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {32,4,8,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {48,4,12,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {64,4,16,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {80,4,20,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {96,4,24,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {112,4,28,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {128,4,32,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {144,4,36,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {160,4,40,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {176,4,44,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {192,8,24,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {224,8,28,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {256,16,16,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {288,8,36,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {320,16,20,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {352,8,44,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {384,16,24,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {448,16,28,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {512,32,16,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {576,16,36,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {640,32,20,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {704,16,44,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {768,32,24,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {896,32,28,0,SmithWatermanKernelConfig::Approach::Unused}, 
            {1024,32,32,0,SmithWatermanKernelConfig::Approach::Unused}, 
        };

        return configs;
    }

    //H100 SXM
    __inline__
    std::vector<SmithWatermanKernelConfig> getOptimalKernelConfigs_SW_sm90(){
        std::vector<SmithWatermanKernelConfig> configs{
            {16,4,4,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {32,4,8,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {48,4,12,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {64,4,16,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {80,4,20,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {96,4,24,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {112,4,28,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {128,4,32,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {144,4,36,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {160,4,40,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {176,4,44,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {192,8,24,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {224,8,28,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {256,8,32,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {288,8,36,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {320,8,40,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {352,8,44,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {384,16,24,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {448,16,28,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {512,16,32,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {576,16,36,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {640,16,40,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {704,16,44,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {768,32,24,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {896,32,28,1,SmithWatermanKernelConfig::Approach::Unused}, 
            {1024,32,32,1,SmithWatermanKernelConfig::Approach::Unused}, 
        };

        return configs;
    }

    __inline__
    std::vector<SmithWatermanKernelConfig> getOptimalKernelConfigs_SW_default(){
        return getOptimalKernelConfigs_SW_sm89();
    }

    __inline__
    std::vector<SmithWatermanKernelConfig> getOptimalKernelConfigs_SW(int deviceId){
        int ccMajor = 0;
        int ccMinor = 0;
        cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, deviceId);
        cudaDeviceGetAttribute(&ccMinor, cudaDevAttrComputeCapabilityMinor, deviceId);

        std::vector<SmithWatermanKernelConfig> configs;
        
        if(ccMajor == 7 && ccMinor == 5){
            configs = getOptimalKernelConfigs_SW_sm75();
        }else if(ccMajor == 8 && ccMinor == 0){
            configs = getOptimalKernelConfigs_SW_sm80();
        }else if(ccMajor == 8 && ccMinor == 9){
            configs = getOptimalKernelConfigs_SW_sm89();
        }else if(ccMajor == 9 && ccMinor == 0){
            configs = getOptimalKernelConfigs_SW_sm90();
        }else{
            configs = getOptimalKernelConfigs_SW_default();
        }

        return configs;
    }


}

#endif