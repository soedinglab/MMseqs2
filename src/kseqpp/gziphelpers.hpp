#ifndef GZIP_HELPERS_HPP
#define GZIP_HELPERS_HPP

#include <zlib.h>
#include <cassert>
#include <fstream>
#include <string>


namespace kseqpp{

    inline
    bool hasGzipHeader(const std::string& filename){
        std::ifstream is(filename, std::ios_base::binary);
        if(!bool(is)){
            throw std::runtime_error("Cannot open file " + filename);
        }
        unsigned char buf[2];
        is.read(reinterpret_cast<char*>(&buf[0]), 2);

        if(buf[0] == 0x1f && buf[1] == 0x8b){
            return true;
        }else{
            return false;
        }
    }

    struct DecompressResult{
        bool mustContinue = false;
        int writtenBytes = 0;
        int statuscode = 0;
    };

    inline
    DecompressResult decompressInMemoryCore(z_stream& zs, unsigned char* output, int outputsize){
        DecompressResult result;

        zs.avail_out = outputsize;
        zs.next_out = output;
        result.statuscode = inflate(&zs, Z_NO_FLUSH);
        assert(result.statuscode != Z_STREAM_ERROR);
        if(result.statuscode < 0){
            return result;
        }

        result.writtenBytes = outputsize - zs.avail_out;
        result.mustContinue = (zs.avail_out == 0);

        return result;
    }

    inline
    DecompressResult decompressInMemory(z_stream& zs, unsigned char* input, int inputsize, unsigned char* output, int outputsize){
        zs.avail_in = inputsize;
        zs.next_in = input;
        return decompressInMemoryCore(zs, output, outputsize);
    }

    inline
    DecompressResult continueDecompressInMemory(z_stream& zs, unsigned char* output, int outputsize){
        return decompressInMemoryCore(zs, output, outputsize);
    }

    inline
    DecompressResult decompressInMemory(z_stream& zs, char* input, int inputsize, char* output, int outputsize){
        return decompressInMemory(zs, 
                                  reinterpret_cast<unsigned char*>(input), inputsize, 
                                  reinterpret_cast<unsigned char*>(output), outputsize);
    }

    inline
    DecompressResult continueDecompressInMemory(z_stream& zs, char* output, int outputsize){
        return continueDecompressInMemory(zs, 
                                          reinterpret_cast<unsigned char*>(output), outputsize);
    }

} // namespace kseqpp

#endif