//
//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#include "ReduceMatrix.h"
#include "SubstitutionMatrix.h"

int main (int argc, const char * argv[])
{
    
    SubstitutionMatrix subMat("../../data/blosum30.out",20);
    
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("%c\t",subMat.int2aa[i]);
    printf("\n");
    ReduceMatrix redMat(subMat.probMatrix,
                        subMat.aa2int,subMat.int2aa,subMat.ALPHABET_SIZE,subMat.ALPHABET_SIZE-16);
    printf("Normal : ");
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("%c\t",subMat.int2aa[i]);
    printf("\nReduced: ");
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("%c\t",redMat.reduced_int2aa[i]);
    printf("\nNormal : ");
    for(int i = 65; i<'Z';i++)
        printf("%d\t",subMat.aa2int[i]); 
    printf("\nReduced: ");
    for(int i = 65; i<'Z';i++)
        printf("%d\t",redMat.reduced_aa2int[i]);    
    printf("\n");

    return 0;
}

