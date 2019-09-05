#include <iostream>
#include <list>
#include <algorithm>
#include <math.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>

#include "Clustering.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Parameters.h"

const char* binary_name = "test_dbreader_zlib";

//static off_t fsize_orDie(const char *filename)
//{
//    struct stat st;
//    if (stat(filename, &st) == 0) return st.st_size;
//    /* error */
//    perror(filename);
//    exit(1);
//}
//
//static FILE* fopen_orDie(const char *filename, const char *instruction)
//{
//    FILE* const inFile = fopen(filename, instruction);
//    if (inFile) return inFile;
//    /* error */
//    perror(filename);
//    exit(2);
//}
//
//static void* malloc_orDie(size_t size)
//{
//    void* const buff = malloc(size);
//    if (buff) return buff;
//    /* error */
//    perror(NULL);
//    exit(3);
//}

//static void* loadFile_orDie(const char* fileName, size_t* size)
//{
//    off_t const fileSize = fsize_orDie(fileName);
//    size_t const buffSize = (size_t)fileSize;
//    if ((off_t)buffSize < fileSize) {   /* narrowcast overflow */
//        fprintf(stderr, "%s : filesize too large \n", fileName);
//        exit(4);
//    }
//    FILE* const inFile = fopen_orDie(fileName, "rb");
//    void* const buffer = malloc_orDie(buffSize);
//    size_t const readSize = fread(buffer, 1, buffSize, inFile);
//    if (readSize != (size_t)buffSize) {
//        fprintf(stderr, "fread: %s : %s \n", fileName, strerror(errno));
//        exit(5);
//    }
//    fclose(inFile);  /* can't fail, read only */
//    *size = buffSize;
//    return buffer;
//}
//
//
//static void saveFile_orDie(const char* fileName, const void* buff, size_t buffSize)
//{
//    FILE* const oFile = fopen_orDie(fileName, "wb");
//    size_t const wSize = fwrite(buff, 1, buffSize, oFile);
//    if (wSize != (size_t)buffSize) {
//        fprintf(stderr, "fwrite: %s : %s \n", fileName, strerror(errno));
//        exit(6);
//    }
//    if (fclose(oFile)) {
//        perror(fileName);
//        exit(7);
//    }
//}


//static void compress_orDie(const char* fname, const char* oname) {
//    size_t fSize;
//    void* const fBuff = loadFile_orDie(fname, &fSize);
//    size_t const cBuffSize = ZSTD_compressBound(fSize);
//    void* const cBuff = malloc_orDie(cBuffSize);
//
//    size_t const cSize = ZSTD_compress(cBuff, cBuffSize, fBuff, fSize, 1);
//    if (ZSTD_isError(cSize)) {
//        fprintf(stderr, "error compressing %s : %s \n", fname, ZSTD_getErrorName(cSize));
//        exit(8);
//    }
//
//    saveFile_orDie(oname, cBuff, cSize);
//
//    /* success */
//    printf("%25s : %6u -> %7u - %s \n", fname, (unsigned)fSize, (unsigned)cSize, oname);
//
//    free(fBuff);
//    free(cBuff);
//}
//
//
//static char* createOutFilename_orDie(const char* filename) {
//    size_t const inL = strlen(filename);
//    size_t const outL = inL + 5;
//    void* const outSpace = malloc_orDie(outL);
//    memset(outSpace, 0, outL);
//    strcat((char*)outSpace, filename);
//    strcat((char*)outSpace, ".zst");
//    return (char*)outSpace;
//}

int main (int, const char**) {
    DBWriter writer("dataLinear", "dataLinear.index", 1, Parameters::WRITER_COMPRESSED_MODE, Parameters::DBTYPE_NUCLEOTIDES);
    writer.open();
    const char * data = "CTGGCGAAACCCAGACCGGTAAGCTTTTCCGTATGCGCGGTAAAGGCGTCAAGTCTGTCC"
                        "GCGGTGGCGCACAGGGTGATTTGCTGTGCCGCGTTGTCGTCGAAACACCGGTAGGCCTGA"
                        "ACGAGAAGCAGAAACAGCTGCTGCAAGAGCTGCAAGAAAGCTTCGGTGGCCCAACCGGTG"
                        "AGCACAACAGCCCGCGCTCAAAGAGCTTCTTTGATGGTGTGAAGAAGTTTTTTGACGACC"
                        "TGACCCGCTAACCTCCCCAAAAGCCTGCCCGTGGGCAGGCCTGGGTAAAAATAGGGTGCG"
                        "TTGAAGATATGCGAGCACCTGTAAAGTGGCGGGGATCACTCCCCGCCGTTGCTCTTACTC"
                        "GGATTCGTAAGCCGTGAAAACAGCAACCTCCGTCTGGCCAGTTCGGGTGTGAACCTCACA"
                        "GAGGTCTTTTCTCGTTACCAGCGCCGCCACTACGGCGGTGATACAGATGACGATCAGGGC"
                        "GACAATCATCGCCTTATGCTGCTTCATTGCTCTCTTCTCCTTGACCTTACGGTCAGTAAG"
                        "AGGCACTCTACATGTGTTCAGCATATAGGGGGCCTCGGGTTGATGGTAAAATATCACTCG"
                        "GGGCTTTTCTCTATCTGCCGTTCAGCTAATGCCTGAGACAGACAGCCTCAAGCACCCGCC"
                        "GCTATTATATCGCTCTCTTTAACCCATTCTGTTTTATCGATTCTAATCCTGAAGACGCCT"
                        "CGCATTTTTATGGCGTAATTTTTTAATGATTTAATTATTTAACTTTAATTTATCTCTTCA";

    writer.writeData((char*)data,strlen(data), 1,0);
    writer.close();
    DBReader<unsigned int> reader("dataLinear", "dataLinear.index", 1, 0);
    reader.open(0);
    reader.readMmapedDataInMemory();
    reader.printMagicNumber();
    std::cout << reader.getSize() << std::endl;
    for(size_t i = 0; i < reader.getSize(); i++){
        std::cout << reader.getSeqLen(i) << std::endl;
        std::cout << reader.getData(i, 0) << std::endl;
    }
    reader.close();

}
