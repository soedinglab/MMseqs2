// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>

#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBConcat.h"
#include "DBWriter.h"
#include "Log.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

enum {
    MSA = 0,
    PSSM,
	ca3m
};

size_t findMaxSetSize(DBReader<unsigned int> *reader) {
    //Find the max set size
    size_t maxSetSize = 0;
    for (size_t i = 0; i < reader->getSize(); i++) {
        char *data = reader->getData(i);
        size_t entry = 0;
        size_t position = 0;
        while (data[position] != '\0') {
            if (data[position] == '\n') {
                entry++;
            }
            position++;
        }
        maxSetSize = std::max(maxSetSize, entry);
    }
    return maxSetSize + 1;
}

MultipleAlignment::MSAResult computeAlignment(MultipleAlignment &aligner, Sequence *centerSequence,
                                              std::vector<Sequence *> seqSet,
                                              std::vector<Matcher::result_t> alnResults,
							      bool allowDeletion, bool sameDatabase) {

	if (alnResults.size() > 0) {
		std::vector<Matcher::result_t> alnWithoutIdentity;
		for (size_t i = 0; i < alnResults.size(); i++) {
			if (alnResults[i].dbKey == centerSequence->getDbKey() && sameDatabase == true) { ;
			} else {
				alnWithoutIdentity.push_back(alnResults[i]);
			}
		}

		return aligner.computeMSA(centerSequence, seqSet, alnWithoutIdentity, !allowDeletion);
	} else {
		return aligner.computeMSA(centerSequence, seqSet, !allowDeletion);
	}
}






int result2outputmode(Parameters &par, int mode) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> *qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
    qDbr->open(DBReader<unsigned int>::NOSORT);

    std::string headerNameQuery(par.db1);
    headerNameQuery.append("_h");

    std::string headerIndexNameQuery(par.db1);
    headerIndexNameQuery.append("_h.index");
	
	std::string headerNameTarget(par.db2);
	headerNameTarget.append("_h");

	std::string headerIndexNameTarget(par.db2);
	headerIndexNameTarget.append("_h.index");

    DBReader<unsigned int> *queryHeaderReader = new DBReader<unsigned int>(headerNameQuery.c_str(), headerIndexNameQuery.c_str());
    queryHeaderReader->open(DBReader<unsigned int>::NOSORT); // NOSORT because the index should be in the same order as resultReader

    DBReader<unsigned int> *tDbr = qDbr;
    DBReader<unsigned int> *tempateHeaderReader = queryHeaderReader;

    size_t maxSequenceLength = 0;
    unsigned int *lengths = qDbr->getSeqLens();
    for (size_t i = 0; i < qDbr->getSize(); i++) {
        maxSequenceLength = std::max((size_t) lengths[i], maxSequenceLength);
    }

    bool sameDatabase = true;
    if (par.db1.compare(par.db2) != 0) {
        sameDatabase = false;
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        tDbr->open(DBReader<unsigned int>::SORT_BY_LINE);

        lengths = tDbr->getSeqLens();
        for (size_t i = 0; i < tDbr->getSize(); i++) {
            maxSequenceLength = std::max((size_t) lengths[i], maxSequenceLength);
        }

        tempateHeaderReader = new DBReader<unsigned int>(headerNameTarget.c_str(), headerIndexNameTarget.c_str());
        tempateHeaderReader->open(DBReader<unsigned int>::NOSORT);
    }
	

	DBConcat *referenceDBr,*referenceHeadersDBr;
	std::string referenceName(par.db4);
	std::string referenceIndexName(par.db4);
	if(mode==ca3m)
	{
		referenceName.append("_ca3m.ffdata");
		referenceIndexName.append("_ca3m.ffindex");
		
		std::string referenceSeqName(par.db4);
		std::string referenceSeqIndexName(par.db4);
		referenceSeqName.append("_sequence.ffdata");
		referenceSeqIndexName.append("_sequence.ffindex");
		
		referenceDBr = new DBConcat(par.db1.c_str(), par.db1Index.c_str(),par.db2.c_str(), par.db2Index.c_str(),referenceSeqName.c_str(),referenceSeqIndexName.c_str(),par.threads);
		referenceDBr->concat();
		referenceDBr->open(DBReader<unsigned int>::SORT_BY_LINE); // When exporting in ca3m,
																		// we need to have an access in SORT_BY_LINE
																		// mode in order to keep track of the original 
																		// line number in the index file.
		
		
		std::string referenceHeadersName(par.db4);
		std::string referenceHeadersIndexName(par.db4);
					
		referenceHeadersName.append("_header.ffdata");
		referenceHeadersIndexName.append("_header.ffindex");
		
		referenceHeadersDBr = new DBConcat(headerNameQuery.c_str(),headerIndexNameQuery.c_str(),headerNameTarget.c_str(), headerIndexNameTarget.c_str(),referenceHeadersName.c_str(),referenceHeadersIndexName.c_str(),par.threads);
		referenceHeadersDBr->concat();
	} else {
			referenceIndexName.append(".index");
	}

    DBReader<unsigned int> *resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str());
    resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);




//    FileUtil::errorIfFileExist(par.db4.c_str());
//    FileUtil::errorIfFileExist(par.db4Index.c_str());

    DBWriter profileWriter(referenceName.c_str(), referenceIndexName.c_str(), par.threads, DBWriter::BINARY_MODE);
    profileWriter.open();
    DBWriter * concensusWriter;
    if(mode == PSSM) {
        concensusWriter = new DBWriter(std::string(par.db4 + "_consensus").c_str(),
                                       std::string(par.db4 + "_consensus.index").c_str(),
                                       par.threads, DBWriter::ASCII_MODE);
        concensusWriter->open();
    }
    size_t maxSetSize = findMaxSetSize(resultReader);
    // adjust score of each match state by -0.2 to trim alignment
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);
    Debug(Debug::INFO) << "Start computing " << (!mode ? "MSAs" : "profiles") << ".\n";
#pragma omp parallel
    {
        Matcher matcher(maxSequenceLength, &subMat, tDbr->getAminoAcidDBSize(), tDbr->getSize(),
                        par.compBiasCorrection);

        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &subMat, &matcher);
        PSSMCalculator calculator(&subMat, maxSequenceLength, par.pca, par.pcb);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat);
        int sequenceType = Sequence::AMINO_ACIDS;
        if (par.profile == true) {
            sequenceType = Sequence::HMM_PROFILE;
        }

        Sequence *centerSequence = new Sequence(maxSequenceLength, subMat.aa2int, subMat.int2aa,
                                                sequenceType, 0, false, par.compBiasCorrection);

#pragma omp for schedule(static)
        for (size_t id = 0; id < resultReader->getSize(); id++) {

            Log::printProgress(id);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char *results = resultReader->getData(id);
            char dbKey[255 + 1];

            // check if already aligned results exists
            std::vector<Matcher::result_t> alnResults;


			 // Get the sequence from the queryDB
            unsigned int queryKey = resultReader->getDbKey(id);
            char *seqData = qDbr->getDataByDBKey(queryKey);
            centerSequence->mapSequence(0, queryKey, seqData);

            std::vector<Sequence *> seqSet;
            while (*results != '\0') {
                Util::parseKey(results, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                double evalue = 0.0;
                char *entry[255];
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // its an aln result
                if(columns >= Matcher::ALN_RES_WITH_OUT_BT_COL_CNT){
                    evalue = strtod (entry[3], NULL);
                }
                // just add sequences if eval < thr. and if key is not the same as the query in case of sameDatabase
                if ( evalue <= par.evalProfile && (key != queryKey || sameDatabase == false) ) {
                    if(columns > Matcher::ALN_RES_WITH_OUT_BT_COL_CNT) {
							Matcher::result_t res = Matcher::parseAlignmentRecord(results);
							alnResults.push_back(res);
                    }
                    const size_t edgeId = tDbr->getId(key);
                    char *dbSeqData = tDbr->getData(edgeId);
                    Sequence *edgeSequence = new Sequence(tDbr->getSeqLens(edgeId), subMat.aa2int, subMat.int2aa,
                                                          Sequence::AMINO_ACIDS, 0, false, false);
                    edgeSequence->mapSequence(0, key, dbSeqData);
                    seqSet.push_back(edgeSequence);
                }
                results = Util::skipLine(results);
            }




			if (alnResults.size() !=  seqSet.size()) // Make sure we have the backtrace info for all the sequences
			{
				//Debug(Debug::INFO) << "BT info is missing... will recompute alignments.\n";
				alnResults.clear(); // will force to recompute alignments together with bt information
			}

			
					
					
            MultipleAlignment::MSAResult res = computeAlignment(aligner, centerSequence, seqSet,
                                                                alnResults, par.allowDeletion, sameDatabase);
            //MultipleAlignment::print(res, &subMat);
			 
			 alnResults = res.alignmentResults;

            MsaFilter::MsaFilterResult filterRes = filter.filter((const char **) res.msaSequence, res.setSize,
                                                                 res.centerLength, static_cast<int>(par.cov * 100),
                                                                 static_cast<int>(par.qid * 100), par.qsc,
                                                                 static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff);

            char *data;
            size_t dataSize;
            std::stringstream msa;
            std::string result;
            switch (mode) {
                case MSA: { // TODO : the first sequence in the MSA seems to be overwritten by the query seq
								// need to allow insertion in the centerSequence

                    for (size_t i = 0; i < filterRes.setSize; i++) {
                        unsigned int key;
                        char *data;
                        if (i == 0) {
                            key = queryKey;
                            data = queryHeaderReader->getDataByDBKey(key);
                        } else {
                            key = seqSet[i - 1]->getDbKey();
                            data = tempateHeaderReader->getDataByDBKey(key);
                        }
                        if (par.addInternalId) {
                            msa << "#" << key << "\n";
                        }
                        msa << ">" << data;// Header
                        for (size_t pos = 0; pos < res.centerLength; pos++) {
                            char aa = filterRes.filteredMsaSequence[i][pos];
                            msa << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int) aa] : '-');
                        }
                        msa << "\n";
                    }
                    result = msa.str();

						//std::cout << result << std::endl;

                    data = (char *) result.c_str();
                    dataSize = result.length();
                }
                    break;


					// This outputs a a3m format --- TODO
/*                case HMM: {


					std::cout<< "====== Cluster " << id << " =====" <<std::endl;

					// First, decide for each column of the alignment
					// if it is a match or insertion state
					  int *matchCount = new int[res.centerLength];
					  for (size_t pos = 0; pos < res.centerLength; pos++) matchCount[pos] = 0;

                    for (size_t i = 0; i < filterRes.setSize; i++) {
                        for (size_t pos = 0; pos < res.centerLength; pos++) {
                             char aa = filterRes.filteredMsaSequence[i][pos];
								matchCount[pos] += (aa < MultipleAlignment::NAA) ? 1:0;
                        }

                    }

					// then write the alignment
					for (size_t i = 0; i < filterRes.setSize; i++) {

						  unsigned int key;
                        char *data;

                        if (i == 0) {
                            key = queryKey;
                            data = queryHeaderReader->getDataByDBKey(key);
                        } else {
                            key = seqSet[i - 1]->getDbKey();
                            data = tempateHeaderReader->getDataByDBKey(key);
                        }
                        if (par.addInternalId) {
                            msa << "#" << key << "\n";
                        }

						  msa << ">" << data;
                        for (size_t pos = 0; pos < res.centerLength; pos++) {

								char aa = filterRes.filteredMsaSequence[i][pos];
								aa = ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int) aa] : '-');

								bool insertionState = matchCount[pos]<filterRes.setSize/2;//centerSequence->int_sequence[pos] > MultipleAlignment::NAA; //
								//std::cout << (std::string)(insertionState ? "GAP":"NNN");
								if (insertionState)
								{
									if (aa!='-')
										msa << (char)tolower(aa); // insert state a2m format : lowercase
									// else do not match an insert state -> nothing written (a2m format)

								} else // else it's a match state and nothing special has to be done
								{
									msa << aa;
								}

                        }
                        msa << "\n";
                    }
                    result = msa.str();

					  	std::cout << result << std::endl;


                    data = (char *) result.c_str();
                    dataSize = result.length();
                }
                    break;
*/

                case ca3m: {	// Here the backtrace information should be present in the alnResults[i].backtrace for all i
					std::stringstream ca3mStream;
					std::stringstream centerSeqStr;
					std::vector<Matcher::result_t> filteredAln; // alignment information for the sequences that passed the filtering step

					// Retrieve the master sequence
					for (size_t pos = 0; pos < centerSequence->L; pos++) {
						 centerSeqStr <<  subMat.int2aa[centerSequence->int_sequence[pos]];
					}

					// Put the query sequence (master sequence) first in the alignment
					Matcher::result_t firstSequence;
					firstSequence.dbKey = queryKey;
					firstSequence.qStartPos = 0;
					firstSequence.dbStartPos = 0;
					firstSequence.backtrace = std::string((centerSeqStr.str()).size(),'M'); // only matches
					filteredAln.push_back(firstSequence);
					
					
					
					
					// Filtering out alnResults using the same filter as previously
					// defined in object filterRes (with sequence coverage, etc.).
					for (size_t i=0;i<alnResults.size();i++)
					{
						if (filterRes.keep[i+1]) // +1 : indice 0 corres. to the query sequence 
							filteredAln.push_back(alnResults.at(i));
					}
					

					// Write the query sequence (== center sequence == consensus sequence)
					msa << ">" << queryHeaderReader->getDataByDBKey(queryKey) << centerSeqStr.str() << "\n;";


					
					int totalNumberofBlocks = 0;
					char zero = 0;
					// then write the alignment
					for (size_t i = 0; i < filteredAln.size(); i++) {
						unsigned int nbOfBlocks = 0;
						std::stringstream blocksDescription,firstBlock;

						  Matcher::result_t aln = filteredAln.at(i);


							// detect the blocks
							for (size_t btIndex = 0; btIndex < (aln.backtrace).size();)
							{
								int indelLen = 0;
								int matchLen = 0;
								char inOrDel;
								// seek the next insertion or deletion
								while(btIndex < (aln.backtrace).size() && aln.backtrace.at(btIndex) == 'M' && matchLen < 255)
								{
									btIndex++;
									matchLen++;
								}

								if (btIndex < (aln.backtrace).size()  && aln.backtrace.at(btIndex) != 'M') // end of block because an I or D was found
									inOrDel = aln.backtrace.at(btIndex); // store whether it is I or D

								// seek the next match
								while(btIndex < (aln.backtrace).size() && aln.backtrace.at(btIndex) == inOrDel && indelLen < 255)
								{
									btIndex++;
									indelLen++;
								}
								// deletion must be count negatively
								// and deletion states are insertion states from the sight of the aligned sequence
								if(indelLen && inOrDel == 'I')
									indelLen *= -1;

								blocksDescription.write(reinterpret_cast<const char*>(&matchLen), sizeof(char));
								blocksDescription.write(reinterpret_cast<const char*>(&indelLen), sizeof(char));
								nbOfBlocks++;
							}

							// add the deletion at the begining of the sequence (local alignment)
							unsigned int firstGap = aln.qStartPos;
							while (firstGap) // need to make as much 127-long deletion blocks as needed
							{
								unsigned int gapToWrite = -std::min((unsigned int)(127),firstGap);
								firstBlock.write(reinterpret_cast<const char*>(&zero),sizeof(char)); // no match in this block,
								firstBlock.write(reinterpret_cast<const char*>(&gapToWrite), sizeof(char)); // only gaps
								firstGap -= -gapToWrite;
								nbOfBlocks++;
							}
							

							
							unsigned short int startPos = aln.dbStartPos + 1; // Starts at 1 in ca3m format
							unsigned int targetIdentifier = referenceDBr->getId(i ? referenceDBr->dbBKeyMap(aln.dbKey) : referenceDBr->dbAKeyMap(aln.dbKey));
							msa.write(reinterpret_cast<const char*>(&targetIdentifier), sizeof(unsigned int));
							msa.write(reinterpret_cast<const char*>(&startPos), sizeof(unsigned short int));
							msa.write(reinterpret_cast<const char*>(&nbOfBlocks), sizeof(unsigned short int));
							msa << firstBlock.str() << blocksDescription.str();


						   totalNumberofBlocks += nbOfBlocks;
                    }

					result = msa.str();
					data = (char *) result.c_str();
					dataSize = result.length();
                }
                    break;
                case PSSM: {
//                    std::cout << centerSequence->getDbKey() << " " << res.setSize << std::endl;;
//                    for (size_t i = 0; i < res.setSize; i++) {
//                        for (size_t pos = 0; pos < res.msaSequenceLength; pos++) {
//                            char aa = res.msaSequence[i][pos];
//                            std::cout << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int) aa] : '-');
//                        }
//                        std::cout << std::endl;
//                    }
//
//                    if (par.pruneSequences) {
//                        filter.pruneAlignment((char **) res.msaSequence, res.setSize, res.centerLength);
//                    }
//
//                    std::cout << centerSequence->getDbKey() << " " << res.setSize << " " << filterRes.setSize <<
//                    std::endl;
//                    for (size_t i = 0; i < filterRes.setSize; i++) {
//                        for (size_t pos = 0; pos < res.msaSequenceLength; pos++) {
//                            char aa = filterRes.filteredMsaSequence[i][pos];
//                            std::cout << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int) aa] : '-');
//                        }
//                        std::cout << std::endl;
//                    }

/*                  std::cout << centerSequence->getDbKey() << " " << res.setSize << " " << filterRes.setSize << std::endl;
		    for (size_t i = 0; i < filterRes.setSize; i++) {
                       for(size_t pos = 0; pos < res.msaSequenceLength; pos++){
                              char aa = filterRes.filteredMsaSequence[i][pos];
                              std::cout << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int)aa] : '-');
                       }
                       std::cout << std::endl;
                    }
*/
		            data = (char *) calculator.computePSSMFromMSA(filterRes.setSize, res.centerLength,
                                                                  filterRes.filteredMsaSequence, par.wg);
                    dataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);

                    std::string consensusStr;
                    for(size_t i = 0; i < res.centerLength; i++){
                        int maxScore = INT_MIN;
                        int maxAA = 0;
                        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
                            if(static_cast<int>(data[i * Sequence::PROFILE_AA_SIZE + aa ]) > maxScore){
                                maxScore = data[i * Sequence::PROFILE_AA_SIZE + aa ];
                                maxAA = aa;
                            }
                        }
                        consensusStr.push_back(subMat.int2aa[maxAA]);
                    }
                    consensusStr.push_back('\n');
                    concensusWriter->write(consensusStr.c_str(), consensusStr.length(), SSTR(queryKey).c_str(), thread_idx);

                    for (size_t i = 0; i < dataSize; i++) {
                        // Avoid a null byte result
                        data[i] = data[i] ^ 0x80;
                    }

                    //pssm.printProfile(res.centerLength);
                    //calculator.printPSSM(res.centerLength);
                }
                    break;
                default:
                    EXIT(EXIT_FAILURE);
            }

            profileWriter.write(data, dataSize, SSTR(queryKey).c_str(), thread_idx);

            MultipleAlignment::deleteMSA(&res);
            for (std::vector<Sequence *>::iterator it = seqSet.begin(); it != seqSet.end(); ++it) {
                Sequence *seq = *it;
                delete seq;
            }
        }
        delete centerSequence;
    }

    // cleanup
    profileWriter.close();
    if(mode == PSSM){
        concensusWriter->close();
        delete concensusWriter;
    }
    resultReader->close();
    delete resultReader;

    if (!sameDatabase) {
        tempateHeaderReader->close();
        tDbr->close();
        delete tempateHeaderReader;
        delete tDbr;
    }

    queryHeaderReader->close();
    qDbr->close();
    delete queryHeaderReader;
    delete qDbr;

    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}


int result2profile(int argc, const char **argv) {
    std::string usage("Calculates profiles from a clustering.\n");
    usage.append("USAGE: <queryDB> <targetDB> <resultDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.result2profile, 4);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::WARNING) << "Compute profile.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);
    int retCode =  result2outputmode(par, PSSM);
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for profile calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    return retCode;
}

int result2msa(int argc, const char **argv) {
    std::string usage("Calculates MSAs from a clustering.\n");
    usage.append("USAGE: <queryDB> <targetDB> <resultDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>\n");
    usage.append("\t & Milot Mirdita <milot@mirdita.de>");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.result2msa, 4);


    Debug(Debug::WARNING) << "Compute Msa.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);
    int retCode;
	if (par.compressMSA)
		retCode = result2outputmode(par, ca3m);
	else
		retCode = result2outputmode(par, MSA);
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for msa calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    return retCode;

}



