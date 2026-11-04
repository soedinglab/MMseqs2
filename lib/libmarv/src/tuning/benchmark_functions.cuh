#ifndef BENCHMARK_FUNCTIONS_CUH
#define BENCHMARK_FUNCTIONS_CUH

#include <chrono>
#include "benchmarkdata.hpp"


std::vector<BenchmarkData> benchmark_float_scoreOnly_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_float_scoreOnly_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

std::vector<BenchmarkData> benchmark_int_scoreOnly_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_int_scoreOnly_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

std::vector<BenchmarkData> benchmark_half2_scoreOnly_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_half2_scoreOnly_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

std::vector<BenchmarkData> benchmark_short2_scoreOnly_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_short2_scoreOnly_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

std::vector<BenchmarkData> benchmark_uint8x4_scoreOnly_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_uint8x4_scoreOnly_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);


std::vector<BenchmarkData> benchmark_float_scoreEndpos_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_float_scoreEndpos_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

std::vector<BenchmarkData> benchmark_int_scoreEndpos_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_int_scoreEndpos_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

std::vector<BenchmarkData> benchmark_half2_scoreEndpos_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_half2_scoreEndpos_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

std::vector<BenchmarkData> benchmark_short2_scoreEndpos_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_short2_scoreEndpos_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

std::vector<BenchmarkData> benchmark_uint8x4_scoreEndpos_singleTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults
);

std::vector<BenchmarkData> benchmark_uint8x4_scoreEndpos_multiTile(
    int timingIterations, 
    std::chrono::milliseconds sleepAmount,
    bool checkResults,
    int tiles
);

#endif