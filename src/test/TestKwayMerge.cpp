//
// Created by mad on 10/26/15.
//

#include <iostream>
#include <queue>


#include "SubstitutionMatrix.h"
#include "Sequence.h"

const char* binary_name = "test_kwaymerge";
struct KmerEntry{
    unsigned int seqId;
    short diagonal;
    KmerEntry(){};
    KmerEntry(unsigned int seqId, short diagonal) : seqId(seqId), diagonal(diagonal){}
};

void mergeKmerEntryLists(KmerEntry **entries, size_t * entrySizes, const int fileCnt);

int main (int, const char**) {
    const int fileCnt = 3;
    std::vector<KmerEntry> array1;
    array1.push_back(KmerEntry(1,0));
    array1.push_back(KmerEntry(2,10));
    array1.push_back(KmerEntry(UINT_MAX,0));
    array1.push_back(KmerEntry(4,0));
    array1.push_back(KmerEntry(8,-2));
    array1.push_back(KmerEntry(UINT_MAX,0));

    std::vector<KmerEntry> array2;
    array2.push_back(KmerEntry(1,0));
    array2.push_back(KmerEntry(2,2));
    array2.push_back(KmerEntry(3,9));
    array2.push_back(KmerEntry(UINT_MAX,0));
    array2.push_back(KmerEntry(5,0));
    array2.push_back(KmerEntry(10,42));
    array2.push_back(KmerEntry(UINT_MAX,0));


    std::vector<KmerEntry> array3;
    array3.push_back(KmerEntry(1,0));
    array3.push_back(KmerEntry(4,2));
    array3.push_back(KmerEntry(UINT_MAX,0));
    array3.push_back(KmerEntry(5,0));
    array3.push_back(KmerEntry(8,42));
    array3.push_back(KmerEntry(UINT_MAX,0));

    KmerEntry ** entries= new KmerEntry*[fileCnt];
    entries[0]=array1.data();
    entries[1]=array2.data();
    entries[2]=array3.data();

    size_t * entryCnts= new size_t[fileCnt];
    entryCnts[0]=array1.size();
    entryCnts[1]=array2.size();
    entryCnts[2]=array3.size();


    mergeKmerEntryLists(entries, entryCnts, fileCnt);
    delete [] entries;
    delete [] entryCnts;
    return 0;
}

struct KmerPosition {
    size_t kmer;
    unsigned int id;
    short pos;
    unsigned int file;

    KmerPosition(){}
    KmerPosition(size_t kmer, unsigned int id,short pos, unsigned int file):
            kmer(kmer), id(id), pos(pos), file(file) {}
};


class CompareResultBySeqId {
public:
    bool operator() (KmerPosition & first, KmerPosition & second) {
        //return (first.eval < second.eval);
        if(first.kmer > second.kmer )
            return true;
        if(second.kmer > first.kmer )
            return false;
        if(first.id > second.id )
            return true;
        if(second.id > first.id )
            return false;
        if(first.pos > second.pos )
            return true;
        if(second.pos > first.pos )
            return false;
        return false;
    }
};

typedef std::priority_queue<KmerPosition, std::vector<KmerPosition>, CompareResultBySeqId> KmerPositionQueue;


size_t queueNextEntry(
        KmerPositionQueue &queue, int file,
        size_t offsetPos,
        KmerEntry *entries, size_t entrySize) {
    if(offsetPos + 1 >= entrySize){
        return offsetPos;
    }
    unsigned int repSeqId = entries[offsetPos].seqId;
    size_t pos = 0;
    while(entries[offsetPos + pos].seqId != UINT_MAX){
        queue.push(KmerPosition(repSeqId, entries[offsetPos+pos].seqId,  entries[offsetPos+pos].diagonal, file));
        pos++;
    }
    queue.push(KmerPosition(repSeqId, UINT_MAX, 0l, file));
    pos++;
    return offsetPos+pos;
}


void mergeKmerEntryLists(KmerEntry **entries, size_t * entrySizes, const int fileCnt) {
    KmerPositionQueue queue;
    size_t * offsetPos = new size_t[fileCnt];
    for(int file = 0; file < fileCnt; file++ ){
        offsetPos[file]=queueNextEntry(queue, file, 0, entries[file], entrySizes[file]);
    }
    KmerPosition prevsKmerPos;
    prevsKmerPos.id = UINT_MAX;
    while(queue.empty() == false) {
        KmerPosition res = queue.top();
        queue.pop();
        if(res.id==UINT_MAX){
            offsetPos[res.file] = queueNextEntry(queue, res.file, offsetPos[res.file],
                                                 entries[res.file], entrySizes[res.file]);
        }
        // if its not a duplicate
        if(prevsKmerPos.id != res.id && res.id!=UINT_MAX){
            std::cout << res.kmer << "\t" << res.id <<"\t"<<res.pos << "\t"<< res.file << std::endl;

        }
        prevsKmerPos = res;

    }
    delete [] offsetPos;
}

