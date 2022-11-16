set(taxonomy_header_files
        taxonomy/NcbiTaxonomy.h
        PARENT_SCOPE
        )


set(taxonomy_source_files
        taxonomy/lca.cpp
        taxonomy/addtaxonomy.cpp
        taxonomy/NcbiTaxonomy.cpp
        taxonomy/filtertaxdb.cpp
        taxonomy/filtertaxseqdb.cpp
        taxonomy/aggregatetax.cpp
        taxonomy/createtaxdb.cpp
        taxonomy/createbintaxonomy.cpp
        taxonomy/createbintaxmapping.cpp
        taxonomy/taxonomyreport.cpp
        taxonomy/TaxonomyExpression.h
        PARENT_SCOPE
        )
