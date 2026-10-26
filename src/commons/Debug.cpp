#include "Debug.h"

#include <iostream>


int Debug::debugLevel = Debug::INFO;



void Debug::setDebugLevel (int i) {
    debugLevel = i;
}

Debug::~Debug(){
    if (level <= ERROR && level <= debugLevel){
        std::cout << std::flush;
        if(interactive){
            std::cerr << "\033[" << Color::FG_RED << "m" << buffer << "\033[" << Color::FG_DEFAULT << "m";;
        }else{
            std::cerr << buffer;
        }
        std::cerr << std::flush;
    } else if(level == WARNING && level <= debugLevel){
        if(interactive){
            std::cout << "\033[" << Color::FG_YELLOW << "m" << buffer << "\033[" << Color::FG_DEFAULT << "m";
        }else{
            std::cout << buffer;
        }
        std::cout << std::flush;
    } else if(level > WARNING && level <= debugLevel) {
        std::cout << buffer;
//            std::cout << std::flush;
    }
}
