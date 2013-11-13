#include "DBReader.h"

DBReader::DBReader(const char* dataFileName_, const char* indexFileName_):
    dataFileName(dataFileName_),
    indexFileName(indexFileName_)
{
    dataSize = 0;
}

void DBReader::open(int sort){
    // open ffindex
    dataFile = fopen(dataFileName, "r");
    indexFile = fopen(indexFileName, "r");

    if( indexFile == NULL) { fferror_print(__FILE__, __LINE__, "DBReader", indexFileName);  exit(EXIT_FAILURE); }
    if( dataFile == NULL) { fferror_print(__FILE__, __LINE__, "DBReader", dataFileName);  exit(EXIT_FAILURE); }

    data = ffindex_mmap_data(dataFile, &dataSize);

    index = ffindex_index_parse(indexFile, 0);

    if(index == NULL)
    {
        fferror_print(__FILE__, __LINE__, "ffindex_index_parse", indexFileName);
        exit(EXIT_FAILURE);
    }

    size = index->n_entries;

    seqLens = new unsigned short [size];

    for (size_t i = 0; i < size; i++){
        seqLens[i] = (unsigned short)(ffindex_get_entry_by_index(index, i)->length);
    }

    // sort sequences by length and generate the corresponding id mappings
    id2local = new size_t[size];
    local2id = new size_t[size];
    for (int i = 0; i < size; i++){
        id2local[i] = i;
        local2id[i] = i;
    }
    
    if (sort == 1){
        calcLocalIdMapping();

        // adapt sequence lengths
        for (size_t i = 0; i < size; i++){
            seqLens[i] = (unsigned short)(ffindex_get_entry_by_index(index, local2id[i])->length);
        }
    }
}

void DBReader::close(){
    fclose(dataFile);
    fclose(indexFile);
    delete[] id2local;
    delete[] local2id;
    delete[] seqLens;
}

char* DBReader::getData (size_t id){
    if (id >= size){
        std::cerr << "getData: local id (" << id << ") >= db size (" << size << ")\n";
        exit(1);
    }
    id = local2id[id];
    if (id >= size){
        std::cerr << "getData: global id (" << id << ") >= db size (" << size << ")\n";
        exit(1);
    }
    if (ffindex_get_entry_by_index(index, id)->offset >= dataSize){ 
        std::cerr << "Invalid database read for id=" << id << "\n";
        exit(1);
    }
    return data + (ffindex_get_entry_by_index(index, id)->offset);
}

char* DBReader::getDataByDBKey (char* key){
    return ffindex_get_data_by_name(data, index, key);
}

size_t DBReader::getSize (){
    return size;
}

char* DBReader::getDbKey (size_t id){
    if (id >= size){
        std::cerr << "getDbKey: local id (" << id << ") >= db size (" << size << ")\n";
        exit(1);
    }

    id = local2id[id];
    if (id >= size){
        std::cerr << "getDbKey: global id (" << id << ") >= db size (" << size << ")\n";
        exit(1);
    }
    return &(ffindex_get_entry_by_index(index, id)->name[0]);
}

unsigned short* DBReader::getSeqLens(){
    return seqLens;
}

void DBReader::merge(size_t* ids, size_t iLeft, size_t iRight, size_t iEnd, size_t* workspace)
{
    size_t i0 = iLeft;
    size_t i1 = iRight;
    size_t j;

    for (j = iLeft; j < iEnd; j++)
    {  
        if (i0 < iRight && (i1 >= iEnd || seqLens[ids[i0]] >= seqLens[ids[i1]]))
        {  
            workspace[j] = ids[i0];
            i0 = i0 + 1;
        }
        else
        {  
            workspace[j] = ids[i1];
            i1 = i1 + 1;
        }
    }
}

void DBReader::sort(size_t* ids, size_t* workspace)
{
    size_t* origIds = ids;
    size_t* tmp;
    for (size_t width = 1; width < size; width = 2 * width)
    {  
        for (size_t i = 0; i < size; i = i + 2 * width)
        {  
            merge(ids, i, std::min(i+width, size), std::min(i+2*width, size), workspace);
        }
        // swap pointers
        tmp = ids;
        ids = workspace;
        workspace = tmp;
    }
    // ensure that the sorted array is stored in the original array
    for (size_t i = 0; i < size; i++)
        origIds[i] = ids[i];
}

/* Sort sequences by length and create two mappings id <-> local id.
 * Necessary for the calculation of the variable sequence similarity threshold in QueryScore.
 * Returns the two mappings id->local and local->id.
 */
void DBReader::calcLocalIdMapping(){
   size_t* workspace = new size_t[size];

   // sort the enties by the length of the sequences
   sort(local2id, workspace);

   for (size_t i = 0; i < size; i++){
       id2local[local2id[i]] = i;
   }

   delete[] workspace;
}
