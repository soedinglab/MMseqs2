#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "KSeqWrapper.h"
#include "MsaFilter.h"
#include "MultipleAlignment.h"

MultipleAlignment::MSAResult readMSA(const std::string &a3mPath) {
    KSeqWrapper* kseq = KSeqFactory(a3mPath.c_str());
    MultipleAlignment::MSAResult msa(0, 0, 0, NULL);
    msa.msaSequenceLength = 0;
    msa.centerLength = 0;
    msa.setSize = 0;
    std::vector<std::string> sequences;
    std::string seq;
    Debug::Progress progress;
    while (kseq->ReadEntry()) {
        progress.updateProgress();
        const KSeqWrapper::KSeqEntry &e = kseq->entry;
        seq.reserve(e.sequence.l);
        if (msa.setSize == 0) {
            msa.msaSequenceLength = e.sequence.l;
            msa.centerLength = e.sequence.l;
        }

        for (size_t i = 0; i < e.sequence.l; i++) {
            char c = e.sequence.s[i];
            if (!(c >= 'a' && c <= 'z')) {
                seq.push_back(c);
            }
        }
        sequences.push_back(seq);
        seq.clear();
        msa.setSize++;
    }
    msa.msaSequence = new char*[msa.setSize];
    for (size_t i = 0; i < msa.setSize; i++) {
        msa.msaSequence[i] = &(sequences[i])[0];
    }
    delete kseq;
    return msa;
}

int filtera3m(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::vector<std::string> qid_str_vec = Util::split(par.qid, ",");
    std::vector<int> qid_vec;
    for (size_t qid_idx = 0; qid_idx < qid_str_vec.size(); qid_idx++) {
        float qid_float = strtod(qid_str_vec[qid_idx].c_str(), NULL);
        qid_vec.push_back(static_cast<int>(qid_float*100));
    }
    std::sort(qid_vec.begin(), qid_vec.end());

    MultipleAlignment::MSAResult msa = readMSA(par.db1.c_str());
    FILE* handle = FileUtil::openFileOrDie(par.db2.c_str(), "w", false);

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, -0.2f);
    MsaFilter filter(msa.centerLength + 1, msa.setSize, &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());
    filter.filter(
        msa.setSize,
        msa.centerLength,
        static_cast<int>(par.covMSAThr * 100),
        qid_vec,
        par.qsc,
        static_cast<int>(par.filterMaxSeqId * 100),
        par.Ndiff,
        par.filterMinEnable,
        (const char **) msa.msaSequence,
        false
    );
    bool* kept = new bool[msa.setSize];
    filter.getKept(kept, msa.setSize);

    KSeqWrapper* kseq2 = KSeqFactory(par.db1.c_str());
    size_t i = 0;
    while (kseq2->ReadEntry()) {
        if (kept[i] == false) {
            i++;
            continue;
        }
        fwrite(">", sizeof(char), 1, handle);
        fwrite(kseq2->entry.name.s, sizeof(char), kseq2->entry.name.l, handle);
        if (kseq2->entry.comment.l > 0) {
            fwrite(" ", sizeof(char), 1, handle);
            fwrite(kseq2->entry.comment.s, sizeof(char), kseq2->entry.comment.l, handle);
        }
        fwrite("\n", sizeof(char), 1, handle);
        fwrite(kseq2->entry.sequence.s, sizeof(char), kseq2->entry.sequence.l, handle);
        fwrite("\n", sizeof(char), 1, handle);
        i++;
    }
    if (fclose(handle) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << par.db2 << "\n";
        return EXIT_FAILURE;
    }

    delete kseq2;
    delete[] kept;
    delete[] msa.msaSequence;

    return EXIT_SUCCESS;
}
