
#include <iostream>
#include <vector>
#include <random>
#include <cstddef>
#include <utility>
#include <string>
#include <algorithm>
#include <string>
#include <thread>
#include <chrono>

#include "common.hpp"
#include "benchmarkdata.hpp"

#include "benchmark_functions.cuh"



int main(int argc, char** argv){

    bool checkResults = false;
    int timingIterations = 10;
    int sleep = 100;
    int tiles = 1; //tiles > 1 -> multitile
    bool printBestConfigs = false;
    bool int8 = false;
    bool half2 = false;
    bool short2 = false;
    bool usefloat = false;
    bool int32 = false;
    bool endpos = false;

    for(int x = 1; x < argc; x++){
        std::string argstring = argv[x];
        if(argstring == "--timingIterations"){
            timingIterations = atoi(argv[x+1]);
            x++;
        }
        if(argstring == "--sleep"){
            sleep = atoi(argv[x+1]);
            x++;
        }
        if(argstring == "--checkResults"){
            checkResults = true;
        }
        if(argstring == "--tiles"){
            tiles = atoi(argv[x+1]);
            x++;
        }
        if(argstring == "--printBestConfigs"){
            printBestConfigs = true;
        }
        if(argstring == "--int8"){
            int8 = true;
        }
        if(argstring == "--int32"){
            int32 = true;
        }        
        if(argstring == "--half2"){
            half2 = true;
        }
        if(argstring == "--short2"){
            short2 = true;
        }
        if(argstring == "--float"){
            usefloat = true;
        }
        if(argstring == "--endpos"){
            endpos = true;
        }
    }


    std::cerr << "checkResults " << checkResults << "\n";
    std::cerr << "timingIterations " << timingIterations << "\n";
    std::cerr << "sleep " << sleep << "\n";
    std::cerr << "tiles " << tiles << "\n";
    std::cerr << "printBestConfigs " << printBestConfigs << "\n";
    std::cerr << "int8 " << int8 << "\n";
    std::cerr << "half2 " << half2 << "\n";
    std::cerr << "short2 " << short2 << "\n";
    std::cerr << "float " << usefloat << "\n";
    std::cerr << "int32 " << int32 << "\n";
    std::cerr << "endpos " << endpos << "\n";


    auto sleepAmount = std::chrono::milliseconds(sleep);

    std::vector<BenchmarkData> allResults;

    if(!endpos){

        if(int8){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_uint8x4_scoreOnly_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_uint8x4_scoreOnly_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }

        if(short2){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_short2_scoreOnly_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_short2_scoreOnly_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }

        if(half2){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_half2_scoreOnly_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_half2_scoreOnly_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }

        if(usefloat){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_float_scoreOnly_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_float_scoreOnly_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }

        if(int32){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_int_scoreOnly_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_int_scoreOnly_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }
    }else{
        if(int8){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_uint8x4_scoreEndpos_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_uint8x4_scoreEndpos_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }

        if(short2){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_short2_scoreEndpos_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_short2_scoreEndpos_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }

        if(half2){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_half2_scoreEndpos_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_half2_scoreEndpos_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }

        if(usefloat){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_float_scoreEndpos_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_float_scoreEndpos_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }

        if(int32){
            if(tiles > 1){
                std::vector<BenchmarkData> vec = benchmark_int_scoreEndpos_multiTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults,
                    tiles
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }else{
                std::vector<BenchmarkData> vec = benchmark_int_scoreEndpos_singleTile(
                    timingIterations, 
                    sleepAmount,
                    checkResults
                );
                allResults.insert(allResults.end(), vec.begin(), vec.end());
            }
        }
    }

    printBenchmarkData_full(allResults);
    
    if(printBestConfigs){
        printBenchmarkData_bestPerTile(allResults);
        printBestAsConfig(allResults);
    }
}