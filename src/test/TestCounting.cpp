#include <iostream>
#include <list>
#include <algorithm>
#include <math.h>
#include <CacheFriendlyOperations.h>
#include <map>


int randBetween(size_t start, size_t end){
    return rand() % (end - start) + start;
}
const size_t ENTRYRANGE = 27000000;
const size_t N = 1000000;


void fillNumbers(unsigned int * data, size_t len) {
    size_t startpos = rand() % (ENTRYRANGE - 30);
    size_t j = 0;

    for (size_t i = 0; i < len; i++) {
        if (i % 20 == 0) {
            startpos = rand() % (ENTRYRANGE - 30);
            j = 0;
        }
        data[i] = std::min(ENTRYRANGE, startpos + j++ + (rand() % 20));
    }
    std::sort(data, data + len);
//    for (size_t i = 0; i < len; i++) {
//        std::cout << data[i] << " ";
//    }
//    std::cout << std::endl;
}

int main(int argc, char **argv)
{
    // DBR;ader test
/*    CacheFriendlyOperations diagonalMatcher(3000000,2);

    unsigned int input[] = {1000001,1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000
                                   ,1000000,2000000,3000000,4000000,4000000};
    CounterResult * output1 = new CounterResult[N];
    CounterResult * output2 = new CounterResult[N];

    size_t resSize = diagonalMatcher.countElements(input, sizeof(input) / sizeof(unsigned int), output1);
    std::cout << resSize << std::endl;
    for(size_t i = 0; i < resSize; i++){
        std::cout << output1[i].id << " " << (int) output1[i].count << std::endl;
    }
    std::cout << std::endl;*/
//
////    // try big data
//    unsigned int * number = new unsigned int[N];
//    fillNumbers(number, N);
//    CacheFriendlyOperations counter2(ENTRYRANGE,N/2048);
//    std::map<int, int> cntMap;
//    for(int i = 0; i < N;i++){
//        cntMap[number[i]]++;
//    }
//    size_t pos = 0;
//    typedef std::map<int, int>::iterator it_type;
//    for(it_type iterator = cntMap.begin(); iterator != cntMap.end(); iterator++) {
//        if(iterator->second > 1){
//            for(size_t i = 1; i < iterator->second; i++){
//                output2[pos].id = iterator->first;
//                pos++;
//            }
//        }
//        // Repeat if you also want to iterate through the second map.
//    }
//    std::sort(output2, output2 + pos);
//    size_t countSize = counter2.countElements(number, N, output1);
//    std::cout << countSize << " " << pos << std::endl;
//    std::sort(output1, output1 + countSize);
//
//
//
//    for(size_t i = 0; i < pos; i++){
//        if(output2[i].id != output1[i].id){
//            std::cout << i << ": " << output1[i] << " " << output2[i] << std::endl;
//        }
//    }
//    delete [] output1;
//    delete [] output2;
//    delete [] number;
}
