#include "blosum.hpp"
#include "util.cuh"

#include <cassert>

namespace cudasw4{

    #ifdef __CUDACC__
    __constant__ std::int8_t deviceBlosum[25*25];
    __constant__ int deviceBlosumDim;
    __constant__ int deviceBlosumDimSquared;
    #endif
    
    std::int8_t hostBlosum[25*25];
    int hostBlosumDim;
    int hostBlosumDimSquared;
    
    //set host and device global variables
    
    
    void setProgramWideBlosum(BlosumType blosumType, const std::vector<int>& deviceIds){
        switch(blosumType){
            case BlosumType::BLOSUM45:
                {
                    const auto blosum = BLOSUM45::get1D();
                    const int dim = BLOSUM45::dim;
                    hostBlosumDim = dim;
                    hostBlosumDimSquared = dim * dim;
                    auto it = std::copy(blosum.begin(), blosum.end(), hostBlosum);
                    assert(std::distance(hostBlosum, it) <= 25 * 25);                
                }
                break;
            case BlosumType::BLOSUM50:
                {
                    const auto blosum = BLOSUM50::get1D();
                    const int dim = BLOSUM50::dim;
                    hostBlosumDim = dim;
                    hostBlosumDimSquared = dim * dim;
                    auto it = std::copy(blosum.begin(), blosum.end(), hostBlosum);
                    assert(std::distance(hostBlosum, it) <= 25 * 25);                
                }
                break;
            case BlosumType::BLOSUM62:
                {
                    const auto blosum = BLOSUM62::get1D();
                    const int dim = BLOSUM62::dim;
                    hostBlosumDim = dim;
                    hostBlosumDimSquared = dim * dim;
                    auto it = std::copy(blosum.begin(), blosum.end(), hostBlosum);
                    assert(std::distance(hostBlosum, it) <= 25 * 25);                
                }
                break;
            case BlosumType::BLOSUM80:
                {
                    const auto blosum = BLOSUM80::get1D();
                    const int dim = BLOSUM80::dim;
                    hostBlosumDim = dim;
                    hostBlosumDimSquared = dim * dim;
                    auto it = std::copy(blosum.begin(), blosum.end(), hostBlosum);
                    assert(std::distance(hostBlosum, it) <= 25 * 25);                
                }
                break;
            case BlosumType::BLOSUM45_20:
                {
                    const auto blosum = BLOSUM45_20::get1D();
                    const int dim = BLOSUM45_20::dim;
                    hostBlosumDim = dim;
                    hostBlosumDimSquared = dim * dim;
                    auto it = std::copy(blosum.begin(), blosum.end(), hostBlosum);
                    assert(std::distance(hostBlosum, it) <= 25 * 25);                
                }
                break;
            case BlosumType::BLOSUM50_20:
                {
                    const auto blosum = BLOSUM50_20::get1D();
                    const int dim = BLOSUM50_20::dim;
                    hostBlosumDim = dim;
                    hostBlosumDimSquared = dim * dim;
                    auto it = std::copy(blosum.begin(), blosum.end(), hostBlosum);
                    assert(std::distance(hostBlosum, it) <= 25 * 25);                
                }
                break;
            case BlosumType::BLOSUM62_20:
                {
                    const auto blosum = BLOSUM62_20::get1D();
                    const int dim = BLOSUM62_20::dim;
                    hostBlosumDim = dim;
                    hostBlosumDimSquared = dim * dim;
                    auto it = std::copy(blosum.begin(), blosum.end(), hostBlosum);
                    assert(std::distance(hostBlosum, it) <= 25 * 25);                
                }
                break;
            case BlosumType::BLOSUM80_20:
                {
                    const auto blosum = BLOSUM80_20::get1D();
                    const int dim = BLOSUM80_20::dim;
                    hostBlosumDim = dim;
                    hostBlosumDimSquared = dim * dim;
                    auto it = std::copy(blosum.begin(), blosum.end(), hostBlosum);
                    assert(std::distance(hostBlosum, it) <= 25 * 25);                
                }
                break;
            default:
                assert(false && "unimplemented blosum copy");
                break;
        }
    #ifdef __CUDACC__
        RevertDeviceId rdi{};
    
        int numGpus = deviceIds.size();
    
        for(int gpu = 0; gpu < numGpus; gpu++){
            cudaSetDevice(deviceIds[gpu]); CUERR;
            cudaMemcpyToSymbol(deviceBlosum, &(hostBlosum[0]), sizeof(std::int8_t) * hostBlosumDim * hostBlosumDim); CUERR;
            cudaMemcpyToSymbol(deviceBlosumDim, &hostBlosumDim, sizeof(int)); CUERR;
            cudaMemcpyToSymbol(deviceBlosumDimSquared, &hostBlosumDimSquared, sizeof(int)); CUERR;
        }
    #endif    
    }

} //namespace cudasw4