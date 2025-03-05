#ifndef LENGTH_PARTITIONS_HPP
#define LENGTH_PARTITIONS_HPP

#include "config.hpp"

#include <array>
#include <limits>

namespace cudasw4{

//length k is in partition i if boundaries[i-1] < k <= boundaries[i]

constexpr auto getLengthPartitionBoundaries(){
 
	constexpr int numLengthPartitions = 36;
	std::array<SequenceLengthT, numLengthPartitions> boundaries{
		48,
		64,
		80,
		96,
		112,
		128,
		144,
		160,
		176,
		192,
		208,
		224,
		240,
		256,
		288,
		320,
		352,
		384,
		416,
		448,
		480,
		512,
		576,
		640,
		704,
		768,
		832,
		896,
		960,
		1024,
		1088,
		1152,
		1216,
		1280,
		8000,
		std::numeric_limits<SequenceLengthT>::max()-1
	};


    return boundaries;
}
    

} //namespace cudasw4

#endif