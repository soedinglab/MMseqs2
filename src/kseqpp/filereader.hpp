#ifndef FILEREADER_HPP
#define FILEREADER_HPP

#include "gziphelpers.hpp"

#include <zlib.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <thread>
#include <mutex>
#include <condition_variable>

namespace kseqpp{


class FileReader{
public:    
    virtual ~FileReader() = default;

    int read(char* outputbuffer, int outputsize){
        return readImpl(outputbuffer, outputsize);
    }
private:
    virtual int readImpl(char* outputbuffer, int outputsize) = 0;
};

class RawReader : public FileReader{
public:    
    RawReader(std::string filename_) 
        : filename(std::move(filename_)),
          inputstream(filename){

        if(!bool(inputstream)){
            throw std::runtime_error("Cannot open file " + filename);
        }
    }

private:
    int readImpl(char* outputbuffer, int outputsize) override{
        inputstream.read(outputbuffer, outputsize);
        return inputstream.gcount();
    }

    std::string filename;
    std::ifstream inputstream;
};

class AsyncRawReader : public FileReader{
public:    
    AsyncRawReader(std::string filename_) 
        : 
          filename(std::move(filename_)),
          inputstream(filename),
          inputthread([&](){inputthreadfunc();}),
          isRunning(true){

        if(!bool(inputstream)){
            throw std::runtime_error("Cannot open file " + filename);
        }
    }

    ~AsyncRawReader(){
        cancel();
        inputthread.join();            
    }

    AsyncRawReader(const AsyncRawReader&) = delete;
    AsyncRawReader(AsyncRawReader&&) = delete;
    AsyncRawReader& operator=(const AsyncRawReader&) = delete;
    AsyncRawReader& operator=(AsyncRawReader&&) = delete;

private:

    struct Buffer{
        static constexpr int bufferSize = 1024 * 64;
        Buffer() : data(std::vector<char>(bufferSize)){}

        int numBytes = 0;
        int begin = 0;
        std::vector<char> data;        
    };

    void cancel(){
        std::unique_lock<std::mutex> ul(commMutex);
        canContinue = false;
        cv_producer.notify_one();
    }

    void inputthreadfunc(){
        Buffer tempbuffer;

        while(canContinue && bool(inputstream)){
            inputstream.read(tempbuffer.data.data(), tempbuffer.data.size());
            tempbuffer.numBytes = inputstream.gcount();
            tempbuffer.begin = 0;

            std::unique_lock<std::mutex> ul(commMutex);
            //std::cerr << "tfunc consumerNeedsNext " << consumerNeedsNext << "\n";
            if(!consumerNeedsNext){
                //std::cerr << "cv_producer.wait\n";
                cv_producer.wait(ul, [&](){return consumerNeedsNext || !canContinue;});
            }
            if(!consumerNeedsNext){
                assert(!canContinue);
                std::cerr << "!canContinue break\n";
                break;
            }
            std::swap(tempbuffer, buffer);
            nextBufferIsReady = true;
            consumerNeedsNext = false;
            //std::cerr << "cv_consumer.notify_one\n";
            cv_consumer.notify_one();
        }

        //std::cerr << canContinue << " " << bool(inputstream) << "\n";

        std::unique_lock<std::mutex> ul(commMutex);
        isRunning = false;
        cv_consumer.notify_one();
    }

    int readImpl(char* outputbuffer, int outputsize) override{
        const int oldOutputsize = outputsize;

        while(outputsize > 0 && (buffer.numBytes > 0 || isRunning)){
            if(buffer.numBytes == 0 && isRunning){
                std::unique_lock<std::mutex> ul(commMutex);
                consumerNeedsNext = true;
                //std::cerr << "set consumerNeedsNext = true\n";
                //std::cerr << "cv_producer.notify_one\n";
                cv_producer.notify_one();

                if(!nextBufferIsReady){                    
                    //std::cerr << "cv_consumer.wait\n";
                    cv_consumer.wait(ul, [&](){return nextBufferIsReady || !isRunning;});
                }
                if(!nextBufferIsReady){
                    std::cerr << "!nextBufferIsReady\n";
                    assert(!isRunning);
                    break;
                }
                nextBufferIsReady = false;

                assert(buffer.begin == 0);
                //std::string s(buffer.data.begin() + buffer.begin, buffer.data.end());
                //std::cout << s ;
            }

            const int bytesToCopy = std::min(outputsize, buffer.numBytes);
            assert(buffer.begin + bytesToCopy <= int(buffer.data.size()));
            std::copy_n(buffer.data.data() + buffer.begin, 
                        bytesToCopy, 
                        outputbuffer + oldOutputsize - outputsize);
            outputsize -= bytesToCopy;
            buffer.numBytes -= bytesToCopy;
            buffer.begin += bytesToCopy;
        }

        if(isRunning){
            assert(oldOutputsize - outputsize > 0);
        }

        return oldOutputsize - outputsize;
    }

    std::string filename;
    std::ifstream inputstream;
    std::thread inputthread;
    bool isRunning = false;
    bool canContinue = true;
    bool nextBufferIsReady = false;
    bool consumerNeedsNext = false;
    Buffer buffer;
    std::mutex commMutex;
    std::condition_variable cv_producer;
    std::condition_variable cv_consumer;
};

class ZlibReader : public FileReader{    
public:    
    ZlibReader(std::string filename_) 
        : filename(std::move(filename_)){

        fh = gzopen(filename.c_str(), "r");

        if(fh == NULL){
            throw std::runtime_error("Cannot open file " + filename);
        }

    }

    ~ZlibReader(){
        gzclose(fh);
    }
private:
    int readImpl(char* outputbuffer, int outputsize) override{
        return gzread(fh, outputbuffer, outputsize);
    }

    std::string filename;
    gzFile fh;
};

template<class RawReader_t>
class GzReaderBase : public FileReader{
public:    
    GzReaderBase(std::string filename_) 
        : 
          filename(std::move(filename_)),
          rawReader(filename),
          compressedBuffer(Buffer(compressedBufferSize)),
          decompressedBuffer(Buffer(decompressedBufferSize)){

        zstream.zalloc = Z_NULL;
        zstream.zfree = Z_NULL;
        zstream.opaque = Z_NULL;
        zstream.avail_in = 0;
        zstream.next_in = Z_NULL;
        int initstatus = inflateInit2(&zstream, 16+MAX_WBITS);
        if (initstatus != Z_OK){
            std::cerr << "Error, inflateInit2 returned " << initstatus << '\n';
        }
        assert(initstatus == Z_OK);
    }

    GzReaderBase(const GzReaderBase&) = delete;
    GzReaderBase(GzReaderBase&&) = delete;
    GzReaderBase& operator=(const GzReaderBase&) = delete;
    GzReaderBase& operator=(GzReaderBase&&) = delete;

private:

    static constexpr int compressedBufferSize = 1024 * 64;
    static constexpr int decompressedBufferSize = 6 * compressedBufferSize;

    struct Buffer{
        Buffer() = default;
        Buffer(int size){
            data.resize(size);
        }

        int numBytes = 0;
        int begin = 0;
        std::vector<char> data;        
    };

    int readImpl(char* outputbuffer, int outputsize) override{
        const int oldOutputsize = outputsize;

        auto iszstreamerror = [&](){
            return decompResult.statuscode < 0;
        };

        auto iszstreamend = [&](){
            return decompResult.statuscode == Z_STREAM_END;
        };

        //can serve request if no error occured, if uncompressed data is available, or if there is compressed data to uncompress
        while(outputsize > 0 && !iszstreamerror() && !(iszstreamend() && decompressedBuffer.numBytes == 0)){

            //if no uncompressed data is available, try to decompress some
            if(decompressedBuffer.numBytes == 0 && !iszstreamerror() && !iszstreamend()){
                
                if(decompResult.mustContinue){
                    //current compressedBuffer is not fully processed yet, continue decompressing the data.
                    decompResult = continueDecompressInMemory(zstream, 
                                                        decompressedBuffer.data.data(), 
                                                        decompressedBuffer.data.size());
                }else{
                    //current compressedBuffer is fully processed, read next chunk from file, then start decompression
                    //TIMERSTARTCPU(fileread);
                    compressedBuffer.numBytes = rawReader.read(reinterpret_cast<char*>(compressedBuffer.data.data()), 
                                                               compressedBuffer.data.size());
                    //TIMERSTOPCPU(fileread);
                    //TIMERSTARTCPU(decomp);
                    if(compressedBuffer.numBytes == 0){
                        decompResult.statuscode = Z_STREAM_END;
                        decompResult.writtenBytes = 0;
                    }else{
                        decompResult = decompressInMemory(zstream, compressedBuffer.data.data(), 
                                                            compressedBuffer.numBytes, 
                                                            decompressedBuffer.data.data(), 
                                                            decompressedBuffer.data.size());
                    }
                    //TIMERSTOPCPU(decomp);
                }
                if(decompResult.statuscode < 0){
                    return decompResult.statuscode;
                }
                decompressedBuffer.numBytes = decompResult.writtenBytes;
                decompressedBuffer.begin = 0;
            }

            const int bytesToCopy = std::min(outputsize, decompressedBuffer.numBytes);
            std::copy_n(decompressedBuffer.data.data() + decompressedBuffer.begin, 
                        bytesToCopy, 
                        outputbuffer + oldOutputsize - outputsize);
            outputsize -= bytesToCopy;
            decompressedBuffer.numBytes -= bytesToCopy;
            decompressedBuffer.begin += bytesToCopy;
        }

        return oldOutputsize - outputsize;
    }

    std::string filename;
    RawReader_t rawReader;
    z_stream zstream; 

    Buffer compressedBuffer;
    Buffer decompressedBuffer;
    DecompressResult decompResult;
};


using GzReader = GzReaderBase<RawReader>;
using AsyncGzReader = GzReaderBase<AsyncRawReader>;


} // namespace kseqpp



#endif