#include "NcbiTaxonomy.h"
#include "Debug.h"

const char* binary_name = "test_taxonomy";

int main (int, const char**) {
    NcbiTaxonomy t("/Users/mirdita/tmp/taxonomy/names.dmp",
                   "/Users/mirdita/tmp/taxonomy/nodes.dmp",
                   "/Users/mirdita/tmp/taxonomy/merged.dmp");
                   //"/Users/mirdita/tmp/taxonomy/delnodes.dmp");
    std::vector<int> taxa;
    taxa.push_back(9);
    taxa.push_back(7);
    TaxonNode const * node = t.LCA(taxa);
    Debug(Debug::INFO) << node->name << "\n";
}
