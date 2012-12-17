//
//
//  Created by Martin Steinegger
//  Copyright (c) 2012 -. All rights reserved.
//
#include "Sequence.h"
#include <iostream>
#include <string>
/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
Sequence::Sequence(size_t maxLen,int* aa2int,char* int2aa)
{
    this->int_sequence = new int[maxLen]; 
    this->aa2int = aa2int;
    this->int2aa = int2aa;
    this->maxLen = maxLen;
    currItPos = -1;
}

/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
Sequence::~Sequence()
{
	delete int_sequence;
}

void Sequence::setId(size_t id){
    this->id = id;
}

void Sequence::mapSequence(const char * sequence){
    this->L=strlen(sequence);
    for(size_t i = 0; i < this->L; i++){
        this->int_sequence[i]=this->aa2int[sequence[i]];
    }
    
}

void Sequence::print() {
    for(int i = 0; i < this->L; i++){
        std::cout << int2aa[this->int_sequence[i]];
    }
    std::cout << std::endl;
    for(int i = 0; i < this->L; i++){
        std::cout << this->int_sequence[i];
    }
    std::cout << std::endl;	 
}

bool Sequence::hasNextKmer(size_t kmerSize) {
    if (((currItPos + 1) + kmerSize) <= this->L)
        return true;
    else
        return false;
}

const int * Sequence::nextKmer(size_t kmerSize) {
    if (hasNextKmer(kmerSize) == true) {
        currItPos += 1;
        return &int_sequence[currItPos];
    } else {
        return 0;
    }
}
