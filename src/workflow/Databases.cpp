#include "Util.h"
#include "Parameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "CommandCaller.h"

#include <cassert>

#include "databases.sh.h"

struct EnvironmentEntry {
    const char* key;
    const char* value;
};

struct DatabaseDownload {
    const char *name;
    const char *description;
    const char *citation;
    const char *url;
    bool hasTaxonomy;
    int dbType;
    const unsigned char *script;
    size_t scriptLength;
    std::vector<EnvironmentEntry> environment;
};

std::vector<DatabaseDownload> downloads = {{
   "UniRef100",
   "The UniProt Reference Clusters provide clustered sets of sequences from the UniProt Knowledgebase.",
   "Suzek et al: UniRef: comprehensive and non-redundant UniProt reference clusters. Bioinformatics 23(10), 1282–1288 (2007)",
   "https://www.uniprot.org/help/uniref",
   true, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
   { }
}, {
    "UniRef90",
    "The UniProt Reference Clusters provide clustered sets of sequences from the UniProt Knowledgebase.",
    "Suzek et al: UniRef: comprehensive and non-redundant UniProt reference clusters. Bioinformatics 23(10), 1282–1288 (2007)",
    "https://www.uniprot.org/help/uniref",
    true, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
    { }
}, {
    "UniRef50",
    "The UniProt Reference Clusters provide clustered sets of sequences from the UniProt Knowledgebase.",
    "Suzek et al: UniRef: comprehensive and non-redundant UniProt reference clusters. Bioinformatics 23(10), 1282–1288 (2007)",
    "https://www.uniprot.org/help/uniref",
    true, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
    { }
}, {
    "UniProtKB",
    "The UniProt Knowledgebase is the central hub for the collection of functional information on proteins, with accurate, consistent and rich annotation.",
    "The UniProt Consortium: UniProt: a worldwide hub of protein knowledge. Nucleic Acids Res 47(D1), D506-515 (2019)",
    "https://www.uniprot.org/help/uniprotkb",
    true, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
    { }
}, {
    "UniProtKB/TrEMBL",
    "UniProtKB/TrEMBL (unreviewed) contains protein sequences associated with computationally generated annotation and large-scale functional characterization.",
    "The UniProt Consortium: UniProt: a worldwide hub of protein knowledge. Nucleic Acids Res 47(D1), D506-515 (2019)",
    "https://www.uniprot.org/help/uniprotkb",
    true, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
    { }
}, {
    "UniProtKB/Swiss-Prot",
    "UniProtKB/Swiss-Prot (reviewed) is a high quality manually annotated and non-redundant protein sequence database, which brings together experimental results, computed features and scientific conclusions.",
    "The UniProt Consortium: UniProt: a worldwide hub of protein knowledge. Nucleic Acids Res 47(D1), D506-515 (2019)",
    "https://uniprot.org",
    true, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
    { }
}, {
    "NR",
    "Non-redundant protein sequences from GenPept, Swissprot, PIR, PDF, PDB, and NCBI RefSeq.",
    "NCBI Resource Coordinators: Database resources of the National Center for Biotechnology Information. Nucleic Acids Res 46(D1), D8-D13 (2018)",
    "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA",
    true, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
    { }
}, {
    "NT",
    "Partially non-redundant nucleotide sequences from all traditional divisions of GenBank, EMBL, and DDBJ excluding GSS, STS, PAT, EST, HTG, and WGS.",
    "NCBI Resource Coordinators: Database resources of the National Center for Biotechnology Information. Nucleic Acids Res 46(D1), D8-D13 (2018)",
    "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA",
    false, Parameters::DBTYPE_NUCLEOTIDES, databases_sh, databases_sh_len,
    { }
}, {
    "GTDB",
    "Genome Taxonomy Database is a phylogenetically consistent, genome-based taxonomy that provides rank-normalized classifications for ~150,000 bacterial and archaeal genomes from domain to genus.",
    "Parks et al: A complete domain-to-species taxonomy for Bacteria and Archaea. Nat Biotechnol 38(9), 1079–1086 (2020)",
    "https://gtdb.ecogenomic.org",
    true, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
    { }
}, {
    "PDB",
    "The Protein Data Bank is the single worldwide archive of structural data of biological macromolecules.",
    "Berman et al: The Protein Data Bank. Nucleic Acids Res 28(1), 235-242 (2000)",
    "https://www.rcsb.org",
    false, Parameters::DBTYPE_AMINO_ACIDS, databases_sh, databases_sh_len,
    { }
}, {
    "PDB70",
    "PDB clustered to 70% sequence identity and enriched using HHblits with Uniclust sequences.",
    "Steinegger et al: HH-suite3 for fast remote homology detection and deep protein annotation. BMC Bioinform 20(1), 473 (2019)",
    "https://github.com/soedinglab/hh-suite",
    false, Parameters::DBTYPE_HMM_PROFILE, databases_sh, databases_sh_len,
    { }
}, {
    "Pfam-A.full",
    "The Pfam database is a large collection of protein families, each represented by multiple sequence alignments and hidden Markov models.",
    "El-Gebali and Mistry et al: The Pfam protein families database in 2019. Nucleic Acids Res 47(D1), D427-D432 (2019)",
    "https://pfam.xfam.org",
    false, Parameters::DBTYPE_HMM_PROFILE, databases_sh, databases_sh_len,
    { }
}, {
    "Pfam-A.seed",
    "The Pfam database is a large collection of protein families, each represented by multiple sequence alignments and hidden Markov models.",
    "El-Gebali and Mistry et al: The Pfam protein families database in 2019. Nucleic Acids Res 47(D1), D427-D432 (2019)",
    "https://pfam.xfam.org",
    false, Parameters::DBTYPE_HMM_PROFILE, databases_sh, databases_sh_len,
    { }
}, {
    "Pfam-B",
    "Pfam-B is a large automatically generated supplement to the Pfam database.",
    "Sonnhammer et al: A new Pfam-B is released. Xfam Blog (2020)",
    "https://xfam.wordpress.com/2020/06/30/a-new-pfam-b-is-released",
    false, Parameters::DBTYPE_HMM_PROFILE, databases_sh, databases_sh_len,
    { }
}, {
    "CDD",
    "Conserved Domain Database is a protein annotation resource consisting of well-annotated MSAs for ancient domains and full-length proteins.",
    "Lu et al: CDD/SPARCLE: the conserved domain database in 2020. Nucleic Acids Res 48(D1), D265–D268 (2020)",
    "https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml",
    false, Parameters::DBTYPE_HMM_PROFILE, databases_sh, databases_sh_len,
    { }
}, {
    "eggNOG",
    "eggNOG is a hierarchical, functionally and phylogenetically annotated orthology resource",
    "Huerta-Cepas et al: eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Res 47(D1), D309–D314 (2019)",
    "http://eggnog5.embl.de",
    false, Parameters::DBTYPE_HMM_PROFILE, databases_sh, databases_sh_len,
    { }
}, {
    "dbCAN2",
    "dbCAN2 is a database of carbohydrate-active enzymes.",
    "Zhang et al: dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. Nucleic Acids Res 46(W1), W95-W101 (2018)",
    "http://bcb.unl.edu/dbCAN2",
    false, Parameters::DBTYPE_HMM_PROFILE, databases_sh, databases_sh_len,
    { }
}, {
    "SILVA",
    "SILVA provides datasets of aligned small and large subunit ribosomal RNA sequences for all three domains of life.",
    "Yilmaz et al: The SILVA and \"All-species Living Tree Project (LTP)\" taxonomic frameworks. Nucleic Acids Res 42(D1), D643-D648 (2014)",
    "https://www.arb-silva.de",
    true, Parameters::DBTYPE_NUCLEOTIDES, databases_sh, databases_sh_len,
    { { "SILVA_REL", "138" } }
}, {
    "Resfinder",
    "ResFinder is a database that captures antimicrobial resistance genes from whole-genome data sets.",
    "Zankari et al: Identification of acquired antimicrobial resistance genes. J Antimicrob Chemother 67(11), 2640-2644 (2012)",
    "https://cge.cbs.dtu.dk/services/ResFinder",
    false, Parameters::DBTYPE_NUCLEOTIDES, databases_sh, databases_sh_len,
    { }
}, {
    "Kalamari",
    "Kalamari contains over 250 genomes chosen to be representative of agents tracked by genome-based foodborne disease surveillance, common contaminants, and diverse phyla and bacterial genera.",
    "Katz et al: Kraken with Kalamari: Contamination Detection. ASM Poster, 270 (2018)",
    "https://github.com/lskatz/Kalamari",
    true, Parameters::DBTYPE_NUCLEOTIDES, databases_sh, databases_sh_len,
    { }
},
};

const int PAD_LEFT = 0;
const int PAD_RIGHT = 1;
void appendPadded(std::string& dst, const std::string& value, size_t n, int direction = PAD_LEFT, char padding = ' ') {
    if (n < value.size()) {
        dst.append(value);
        return;
    }
    if (direction == PAD_RIGHT) {
        dst.append(n - value.size(), padding);
    }
    dst.append(value);
    if (direction == PAD_LEFT) {
        dst.append(n - value.size(), padding);
    }
}

std::string listDatabases(const Command &command, bool detailed) {
    size_t nameWidth = 4, urlWidth = 3, dbTypeWidth = 4;
    for (size_t i = 0; i < downloads.size(); ++i) {
        nameWidth = std::max(nameWidth, strlen(downloads[i].name));
        urlWidth = std::max(urlWidth, strlen(downloads[i].url));
        dbTypeWidth = std::max(dbTypeWidth, strlen(Parameters::getDbTypeName(downloads[i].dbType)));
    }

    std::string description;
    description.reserve(1024);
    if (detailed) {
        description += " By ";
        description += command.author;
        description += "\n";
    }

    description += "\n  ";
    appendPadded(description, "Name", nameWidth);
    description.append(1, '\t');
    appendPadded(description, "Type", dbTypeWidth);
    description.append(1, '\t');
    appendPadded(description, "Taxonomy", 8);
    description.append(1, '\t');
    appendPadded(description, "Url", urlWidth);
    description.append(1, '\n');

    for (size_t i = 0; i < downloads.size(); ++i) {
        description.append("- ");
        appendPadded(description, downloads[i].name, nameWidth);
        description.append(1, '\t');
        appendPadded(description, Parameters::getDbTypeName(downloads[i].dbType), dbTypeWidth);
        description.append(1, '\t');
        appendPadded(description, (downloads[i].hasTaxonomy ? "yes" : "-"), 8, PAD_RIGHT);
        description.append(1, '\t');
        // last field in line should not be padded
        //appendPadded(description, downloads[i].url, urlWidth);
        description.append(downloads[i].url);
        description.append(1, '\n');
        if (detailed) {
            if (strlen(downloads[i].description) > 0) {
                description.append(2, ' ');
                description.append(downloads[i].description);
                description.append(1, '\n');
            }
            if (strlen(downloads[i].citation) > 0) {
                description.append("  Cite: ");
                description.append(downloads[i].citation);
                description.append(1, '\n');
            }
            description.append(1, '\n');
        }
    }

    return description;
}

int databases(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    std::string description = listDatabases(command, par.help);
    if (par.filenames.size() == 0 || par.help) {
        par.printUsageMessage(command, par.help ? MMseqsParameter::COMMAND_EXPERT : 0, description.c_str());
        EXIT(EXIT_SUCCESS);
    }

    ssize_t downloadIdx = -1;
    for (size_t i = 0; i < downloads.size(); ++i) {
        if (par.db1 == std::string(downloads[i].name)) {
            downloadIdx = i;
            break;
        }
    }
    if (downloadIdx == -1) {
        par.printUsageMessage(command, par.help ? MMseqsParameter::COMMAND_EXPERT : 0, description.c_str());
        Debug(Debug::ERROR) << "Selected database " << par.db1 << " was not found\n";
        EXIT(EXIT_FAILURE);
    }
    par.printParameters(command.cmd, argc, argv, par.databases);
    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.databases));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    for (size_t i = 0; i < downloads[downloadIdx].environment.size(); ++i) {
        cmd.addVariable(downloads[downloadIdx].environment[i].key, downloads[downloadIdx].environment[i].value);
    }
    cmd.addVariable("TAXONOMY", downloads[downloadIdx].hasTaxonomy ? "TRUE" : NULL);
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERB_PAR", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("COMP_PAR", par.createParameterString(par.verbandcompression).c_str());
    // aria2c gives an (undocumented error with more than 16 connections)
    cmd.addVariable("ARIA_NUM_CONN", SSTR(std::min(16, par.threads)).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    std::string program = tmpDir + "/download.sh";
    FileUtil::writeFile(program, downloads[downloadIdx].script, downloads[downloadIdx].scriptLength);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    EXIT(EXIT_FAILURE);
}
