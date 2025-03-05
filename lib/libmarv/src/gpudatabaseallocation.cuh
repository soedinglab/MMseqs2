#ifndef GPU_DATABASE_ALLOCATION_CUH
#define GPU_DATABASE_ALLOCATION_CUH

#include "config.hpp"

#include "hpc_helpers/simple_allocation.cuh"

struct GpuDatabaseAllocationBase{
    virtual const char* getCharData() const = 0;
    virtual const SequenceLengthT* getLengthData() const = 0;
    virtual const size_t* getOffsetData() const = 0;

    virtual char* getCharData() = 0;
    virtual SequenceLengthT* getLengthData() = 0;
    virtual size_t* getOffsetData() = 0;

    virtual size_t getNumChars() const = 0;
    virtual size_t getNumSubjects() const = 0;
};


struct GpuDatabaseAllocation : public GpuDatabaseAllocationBase{

    GpuDatabaseAllocation() : GpuDatabaseAllocation(0,0) {}
    GpuDatabaseAllocation(size_t numChars, size_t numSubjects){
        d_fulldb_chardata.resize(numChars);
        d_fulldb_lengthdata.resize(numSubjects);
        d_fulldb_offsetdata.resize(numSubjects+1);
    }

    const char* getCharData() const override{
        return d_fulldb_chardata.data();
    }
    const SequenceLengthT* getLengthData() const override{
        return d_fulldb_lengthdata.data();
    }
    const size_t* getOffsetData() const override{
        return d_fulldb_offsetdata.data();
    }
    char* getCharData() override{
        return d_fulldb_chardata.data();
    }
    SequenceLengthT* getLengthData() override{
        return d_fulldb_lengthdata.data();
    }
    size_t* getOffsetData() override{
        return d_fulldb_offsetdata.data();
    }

    size_t getNumChars() const override{
        return d_fulldb_chardata.size();
    }

    size_t getNumSubjects() const override{
        return d_fulldb_lengthdata.size();
    }

    helpers::SimpleAllocationDevice<char, 0> d_fulldb_chardata;
    helpers::SimpleAllocationDevice<SequenceLengthT, 0> d_fulldb_lengthdata;
    helpers::SimpleAllocationDevice<size_t, 0> d_fulldb_offsetdata;
};

struct GpuDatabaseAllocationView : public GpuDatabaseAllocationBase{

    GpuDatabaseAllocationView() = default;
    GpuDatabaseAllocationView(
        char* chardata_, 
        SequenceLengthT* lengthdata_,
        size_t* offsetdata_,
        size_t numChars_,
        size_t numSubjects_
    ): chardata(chardata_), lengthdata(lengthdata_), offsetdata(offsetdata_),numChars(numChars_), numSubjects(numSubjects_){

    }

    const char* getCharData() const override{
        return chardata;
    }
    const SequenceLengthT* getLengthData() const override{
        return lengthdata;
    }
    const size_t* getOffsetData() const override{
        return offsetdata;
    }
    char* getCharData() override{
        return chardata;
    }
    SequenceLengthT* getLengthData() override{
        return lengthdata;
    }
    size_t* getOffsetData() override{
        return offsetdata;
    }

    size_t getNumChars() const override{
        return numChars;
    }

    size_t getNumSubjects() const override{
        return numSubjects;
    }

    char* chardata; //numChars
    SequenceLengthT* lengthdata; // numSubjects
    size_t* offsetdata; //numSubjects + 1

    size_t numChars;
    size_t numSubjects;
};



#endif