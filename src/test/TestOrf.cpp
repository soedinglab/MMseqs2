#include <iostream>
#include <list>
#include <algorithm>
#include <math.h>
#include <SetCover2.h>
#include <Orf.h>


#include "Clustering.h"
#include "SetElement.h"

#include "DBReader.h"
#include "DBWriter.h"



int main(int argc, char **argv)
{

    const char * seq = "CGAAGCGGGTGATGGCCGGCGCCGCGCCGGTTGGCGGCTGGCCATTCAAGGAGTTGAGGAGATGGTCACTGGGCAGCGCGCCGGGGGGCGGCAGCAGCCCAAGGGTCGGGTCATTCCCGATTGGCCGCACCAGGCGCCCGCCACAGCCGGA";
    Orf orf(seq);
    std::vector<Orf::SequenceLocation> res;
    size_t orfMinLength = 1;
    size_t orfMaxLength = SIZE_MAX;
    size_t orfMaxGaps = SIZE_MAX;
    bool orfSkipIncomplete = false;

    orf.FindOrfs(res, orfMinLength, orfMaxLength, orfMaxGaps);

    size_t orfs_num = 0;
    printf("Orf found: %d\n", res.size());

    for (std::vector<Orf::SequenceLocation>::const_iterator it = res.begin(); it != res.end(); it++) {

        Orf::SequenceLocation loc = (Orf::SequenceLocation) * it;

        if (orfSkipIncomplete && (loc.uncertainty_from != Orf::UNCERTAINTY_UNKOWN || loc.uncertainty_to != Orf::UNCERTAINTY_UNKOWN))
            continue;

        printf("%s [Orf: %zu, %zu, %d, %d, %d]\n", "Test", loc.from, loc.to, loc.strand, loc.uncertainty_from, loc.uncertainty_to);

        std::unique_ptr<char[]> sequence(orf.View(loc));
        char* seq = sequence.get();
        size_t length = strlen(seq);
        printf("%s\n", seq);

        orfs_num++;
    }
}
