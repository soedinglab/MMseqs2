#ifndef BLOSUM_HPP
#define BLOSUM_HPP

#include "types.hpp"
#include "util.cuh"
#include <array>
#include <string>
#include <vector>
#include <cstdint>

namespace cudasw4{

#ifdef __CUDACC__

extern __constant__ std::int8_t deviceBlosum[25*25];
extern __constant__ int deviceBlosumDim;
extern __constant__ int deviceBlosumDimSquared;

#endif

extern std::int8_t hostBlosum[25*25];
extern int hostBlosumDim;
extern int hostBlosumDimSquared;

//set host and device global blosum variables
void setProgramWideBlosum(BlosumType blosumType, const std::vector<int>& deviceIds);

} //namespace cudasw4

#endif