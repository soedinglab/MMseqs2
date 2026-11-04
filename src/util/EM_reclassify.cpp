#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "FastSort.h"
#include "Timer.h"
#include "EMTargetGroups.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef OPENMP
#include <omp.h>
#endif

namespace {

const double LOG_COMPATIBILITY_MIN = -60.0;
const double LOG_FRAGMENT_MIN = -120.0;
const double STEP_LENGTH_MIN = 1.0;
const double STEP_LENGTH_CAP = 1.0e4;
const double STEP_LENGTH_GROWTH = 4.0;
const int MIN_ITERATIONS = 6;
const unsigned int NONE = EMTargetGroups::NONE;

int usableThreads(int threads) {
#ifdef OPENMP
    return std::max(1, threads);
#else
    (void)threads;
    return 1;
#endif
}

size_t countRecords(const char *data) {
    size_t records = 0;
    char *cursor = const_cast<char *>(data);
    while (*cursor != '\0') {
        ++records;
        cursor = Util::skipLine(cursor);
    }
    return records;
}

// read name -> fragment name: the mate suffix /1 or /2 is removed, so 'x/1' and 'x/2' (and two reads both named
// 'x', as createdb names Illumina 1.8+ mates) form one fragment
std::string fragmentOfReadName(const std::string &name) {
    const size_t n = name.size();
    if (n > 2 && name[n - 2] == '/' && (name[n - 1] == '1' || name[n - 1] == '2')) {
        return name.substr(0, n - 2);
    }
    return name;
}

struct Model {
    // reads (DB entries with hits) and their hits
    std::vector<unsigned int> queryKey;
    std::vector<size_t> queryOffset;
    std::vector<int> hitScore;
    std::vector<unsigned int> hitGroup;
    std::vector<double> hitPhi;              // phi_qp
    std::vector<double> hitPhiQG;            // sum of phi over the read's hits in the same group
    std::vector<Matcher::result_t> hitResult;
    bool hasBacktrace = false;
    bool hasOrfPosition = false;

    // fragments: queries of a fragment, CSR
    std::vector<size_t> unitOffset;
    std::vector<unsigned int> unitQuery;
    std::vector<unsigned int> unitOfQuery;
    // fragment x group compatibilities, CSR over fragments
    std::vector<size_t> ugOffset;
    std::vector<unsigned int> ugGroup;
    std::vector<double> ugPhi;

    std::vector<double> groupLength;

    size_t units() const { return unitOffset.size() - 1; }
    size_t groupCount() const { return groupLength.size(); }
};

// ------------------------------------------------------------------------------------------- loading
void loadHits(DBReader<unsigned int> &reader, Model &m, int threads) {
    const int nThreads = usableThreads(threads);
    std::vector<size_t> entryIds;
    for (size_t i = 0; i < reader.getSize(); ++i) {
        if (reader.getEntryLen(i) > 1) {
            entryIds.push_back(i);
        }
    }
    std::vector<size_t> records(entryIds.size(), 0);
#pragma omp parallel for schedule(dynamic, 64) num_threads(nThreads)
    for (size_t q = 0; q < entryIds.size(); ++q) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        records[q] = countRecords(reader.getData(entryIds[q], thread_idx));
    }
    m.queryKey.resize(entryIds.size());
    m.queryOffset.assign(entryIds.size() + 1, 0);
    for (size_t q = 0; q < entryIds.size(); ++q) {
        m.queryKey[q] = reader.getDbKey(entryIds[q]);
        m.queryOffset[q + 1] = m.queryOffset[q] + records[q];
    }
    const size_t hitCount = m.queryOffset.back();
    m.hitScore.resize(hitCount);
    m.hitGroup.resize(hitCount);
    m.hitResult.resize(hitCount);
    m.hitPhi.resize(hitCount);
    m.hitPhiQG.resize(hitCount);
    std::vector<char> sawBt(nThreads, 0), sawOrf(nThreads, 0);
#pragma omp parallel num_threads(nThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        const char *entry[255];
#pragma omp for schedule(dynamic, 64)
        for (size_t q = 0; q < entryIds.size(); ++q) {
            char *data = reader.getData(entryIds[q], thread_idx);
            size_t h = m.queryOffset[q];
            while (*data != '\0') {
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    Debug(Debug::ERROR) << "Invalid alignment result record in query " << m.queryKey[q] << ".\n";
                    EXIT(EXIT_FAILURE);
                }
                if (columns == Matcher::ALN_RES_WITH_BT_COL_CNT || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT) {
                    sawBt[thread_idx] = 1;
                }
                if (columns == Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT) {
                    sawOrf[thread_idx] = 1;
                }
                Matcher::result_t r = Matcher::parseAlignmentRecord(data, true);
                m.hitScore[h] = r.score;
                m.hitResult[h] = std::move(r);
                ++h;
                data = Util::skipLine(data);
            }
        }
    }
    for (int t = 0; t < nThreads; ++t) {
        m.hasBacktrace = m.hasBacktrace || sawBt[t];
        m.hasOrfPosition = m.hasOrfPosition || sawOrf[t];
    }
}

// hit -> group.  Genomes (genomeOfKey, --ref-group 1) come first with their summed length over the whole target DB;
// every other hit target is a group of its own, numbered after them, with its own length.  Targets without a hit
// get no group, so a large target DB costs nothing without --ref-group.
void assignGroups(Model &m, const std::vector<unsigned int> &genomeOfKey, const std::vector<double> &genomeLength,
                  int threads) {
    const size_t H = m.hitResult.size();
    unsigned int maxKey = 0;
#pragma omp parallel for schedule(static) reduction(max:maxKey) num_threads(usableThreads(threads))
    for (size_t h = 0; h < H; ++h) {
        maxKey = std::max(maxKey, m.hitResult[h].dbKey);
    }
    m.groupLength = genomeLength;
    std::vector<unsigned int> ownGroup(H > 0 ? static_cast<size_t>(maxKey) + 1 : 0, NONE);
    for (size_t h = 0; h < H; ++h) {
        const unsigned int key = m.hitResult[h].dbKey;
        if (key < genomeOfKey.size() && genomeOfKey[key] != NONE) {
            m.hitGroup[h] = genomeOfKey[key];
            continue;
        }
        if (ownGroup[key] == NONE) {
            ownGroup[key] = static_cast<unsigned int>(m.groupLength.size());
            m.groupLength.push_back(static_cast<double>(std::max(1u, m.hitResult[h].dbLen)));
        }
        m.hitGroup[h] = ownGroup[key];
    }
}

// fragments: every read alone, or (--mate-pair 1) the reads whose .lookup names agree up to a trailing /1 or /2
void buildFragments(Model &m, const std::string &queryDb, bool matePair) {
    const size_t Q = m.queryKey.size();
    m.unitOfQuery.resize(Q);
    for (size_t q = 0; q < Q; ++q) {
        m.unitOfQuery[q] = static_cast<unsigned int>(q);
    }
    size_t U = Q;
    if (matePair) {
        std::unordered_map<unsigned int, unsigned int> qOfKey;
        qOfKey.reserve(Q);
        for (size_t q = 0; q < Q; ++q) {
            qOfKey[m.queryKey[q]] = static_cast<unsigned int>(q);
        }
        std::unordered_map<std::string, unsigned int> fragOfName;
        fragOfName.reserve(Q);
        std::vector<char> named(Q, 0);
        U = 0;
        const bool hasLookup = EMTargetGroups::forEachLookup(queryDb, [&](unsigned int key, const std::string &name, unsigned int) {
            std::unordered_map<unsigned int, unsigned int>::const_iterator q = qOfKey.find(key);
            if (q == qOfKey.end()) {
                return;
            }
            std::pair<std::unordered_map<std::string, unsigned int>::iterator, bool> it =
                    fragOfName.emplace(fragmentOfReadName(name), static_cast<unsigned int>(U));
            U += it.second ? 1 : 0;
            m.unitOfQuery[q->second] = it.first->second;
            named[q->second] = 1;
        });
        if (!hasLookup) {
            Debug(Debug::ERROR) << queryDb << ".lookup is missing: --mate-pair 1 needs the read names.\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t q = 0; q < Q; ++q) {
            if (!named[q]) {
                m.unitOfQuery[q] = static_cast<unsigned int>(U++);     // not in the .lookup: a fragment of its own
            }
        }
    }
    m.unitOffset.assign(U + 1, 0);
    for (size_t q = 0; q < Q; ++q) {
        m.unitOffset[m.unitOfQuery[q] + 1]++;
    }
    for (size_t u = 0; u < U; ++u) {
        m.unitOffset[u + 1] += m.unitOffset[u];
    }
    m.unitQuery.resize(Q);
    std::vector<size_t> fill(m.unitOffset.begin(), m.unitOffset.end() - 1);
    for (size_t q = 0; q < Q; ++q) {
        m.unitQuery[fill[m.unitOfQuery[q]]++] = static_cast<unsigned int>(q);
    }
}

// phi_qp, phi_qg and the fragment x group compatibilities
void buildCompatibilities(Model &m, double lambda, int threads) {
    const size_t Q = m.queryKey.size();
#pragma omp parallel for schedule(dynamic, 256) num_threads(usableThreads(threads))
    for (size_t q = 0; q < Q; ++q) {
        const size_t b = m.queryOffset[q], e = m.queryOffset[q + 1];
        int smax = std::numeric_limits<int>::min();
        for (size_t h = b; h < e; ++h) {
            smax = std::max(smax, m.hitScore[h]);
        }
        for (size_t h = b; h < e; ++h) {
            m.hitPhi[h] = std::exp(std::max(LOG_COMPATIBILITY_MIN, -lambda * static_cast<double>(smax - m.hitScore[h])));
        }
        // phi_qg: hits of one read are few, a quadratic pass is cheaper than a map
        for (size_t h = b; h < e; ++h) {
            double s = 0.0;
            for (size_t k = b; k < e; ++k) {
                if (m.hitGroup[k] == m.hitGroup[h]) {
                    s += m.hitPhi[k];
                }
            }
            m.hitPhiQG[h] = s;
        }
    }

    const size_t U = m.units();
    std::vector<std::vector<std::pair<unsigned int, double> > > perUnit(U);
    size_t keptLenient = 0, withCommon = 0;
#pragma omp parallel num_threads(usableThreads(threads))
    {
        std::vector<std::pair<unsigned int, double> > logPhi;   // (group, log phi_qg) over the unit's reads
        std::vector<std::pair<unsigned int, std::pair<double, unsigned int> > > merged;
#pragma omp for schedule(dynamic, 256) reduction(+:keptLenient, withCommon)
        for (size_t u = 0; u < U; ++u) {
            logPhi.clear();
            unsigned int matesHit = 0;
            for (size_t j = m.unitOffset[u]; j < m.unitOffset[u + 1]; ++j) {
                const unsigned int q = m.unitQuery[j];
                const size_t b = m.queryOffset[q], e = m.queryOffset[q + 1];
                if (b == e) {
                    continue;
                }
                ++matesHit;
                for (size_t h = b; h < e; ++h) {
                    bool first = true;
                    for (size_t k = b; k < h; ++k) {
                        if (m.hitGroup[k] == m.hitGroup[h]) {
                            first = false;
                            break;
                        }
                    }
                    if (first) {
                        logPhi.emplace_back(m.hitGroup[h], std::log(m.hitPhiQG[h]));
                    }
                }
            }
            std::sort(logPhi.begin(), logPhi.end());
            merged.clear();
            for (size_t i = 0; i < logPhi.size(); ++i) {
                if (!merged.empty() && merged.back().first == logPhi[i].first) {
                    merged.back().second.first += logPhi[i].second;
                    merged.back().second.second += 1;
                } else {
                    merged.emplace_back(logPhi[i].first, std::make_pair(logPhi[i].second, 1u));
                }
            }
            bool anyFull = false;
            for (size_t i = 0; i < merged.size(); ++i) {
                anyFull = anyFull || merged[i].second.second == matesHit;
            }
            if (anyFull && matesHit > 1) {
                ++withCommon;
            } else if (matesHit > 1) {
                ++keptLenient;
            }
            double mx = -std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < merged.size(); ++i) {
                if (!anyFull || merged[i].second.second == matesHit) {
                    mx = std::max(mx, merged[i].second.first);
                }
            }
            std::vector<std::pair<unsigned int, double> > &out = perUnit[u];
            for (size_t i = 0; i < merged.size(); ++i) {
                if (anyFull && merged[i].second.second != matesHit) {
                    continue;
                }
                out.emplace_back(merged[i].first, std::exp(std::max(LOG_FRAGMENT_MIN, merged[i].second.first - mx)));
            }
        }
    }
    m.ugOffset.assign(U + 1, 0);
    for (size_t u = 0; u < U; ++u) {
        m.ugOffset[u + 1] = m.ugOffset[u] + perUnit[u].size();
    }
    m.ugGroup.resize(m.ugOffset.back());
    m.ugPhi.resize(m.ugOffset.back());
    for (size_t u = 0; u < U; ++u) {
        for (size_t i = 0; i < perUnit[u].size(); ++i) {
            m.ugGroup[m.ugOffset[u] + i] = perUnit[u][i].first;
            m.ugPhi[m.ugOffset[u] + i] = perUnit[u][i].second;
        }
    }
    Debug(Debug::INFO) << U << " fragments, " << m.ugOffset.back() << " fragment x group candidates";
    if (withCommon + keptLenient > 0) {
        Debug(Debug::INFO) << "; " << withCommon << " mate pairs with a group common to both mates, "
                           << keptLenient << " without (their mates' candidates are kept)";
    }
    Debug(Debug::INFO) << ".\n";
}

// ------------------------------------------------------------------------------------------- EM
struct EmState {
    std::vector<double> theta;
    std::vector<double> gamma;        // per fragment x group candidate
    std::vector<char> active;         // per group: has a candidate
};

// one EM map x -> M(x); also fills gamma
std::vector<double> emMap(const Model &m, EmState &s, const std::vector<double> &x, const std::vector<double> &invLen,
                          double priorTerm, int threads) {
    const size_t U = m.units(), G = m.groupCount();
    std::vector<double> dens(G);
    for (size_t g = 0; g < G; ++g) {
        dens[g] = s.active[g] ? x[g] * invLen[g] : 0.0;
    }
    const int nThreads = usableThreads(threads);
    std::vector<std::vector<double> > partial(nThreads, std::vector<double>(G, 0.0));
#pragma omp parallel num_threads(nThreads)
    {
        int t = 0;
#ifdef OPENMP
        t = omp_get_thread_num();
#endif
        std::vector<double> &acc = partial[t];
#pragma omp for schedule(static)
        for (size_t u = 0; u < U; ++u) {
            double den = 0.0;
            for (size_t i = m.ugOffset[u]; i < m.ugOffset[u + 1]; ++i) {
                const double w = m.ugPhi[i] * dens[m.ugGroup[i]];
                s.gamma[i] = w;
                den += w;
            }
            for (size_t i = m.ugOffset[u]; i < m.ugOffset[u + 1]; ++i) {
                s.gamma[i] = den > 0.0 ? s.gamma[i] / den : 0.0;
                acc[m.ugGroup[i]] += s.gamma[i];
            }
        }
    }
    std::vector<double> next(G, 0.0);
    double sum = 0.0;
    for (size_t g = 0; g < G; ++g) {
        if (!s.active[g]) {
            continue;
        }
        double r = priorTerm;
        for (int t = 0; t < nThreads; ++t) {
            r += partial[t][g];
        }
        next[g] = r;
        sum += r;
    }
    if (sum > 0.0) {
        for (size_t g = 0; g < G; ++g) {
            next[g] /= sum;
        }
    }
    return next;
}

double logLikelihood(const Model &m, const EmState &s, const std::vector<double> &x, const std::vector<double> &invLen,
                     int threads) {
    const size_t U = m.units();
    double ll = 0.0;
    size_t n = 0;
#pragma omp parallel for schedule(static) reduction(+:ll, n) num_threads(usableThreads(threads))
    for (size_t u = 0; u < U; ++u) {
        double mix = 0.0;
        for (size_t i = m.ugOffset[u]; i < m.ugOffset[u + 1]; ++i) {
            const unsigned int g = m.ugGroup[i];
            if (s.active[g]) {
                mix += m.ugPhi[i] * x[g] * invLen[g];
            }
        }
        if (mix > 0.0) {
            ll += std::log(mix);
            ++n;
        }
    }
    return n > 0 ? ll / static_cast<double>(n) : 0.0;
}

void projectSimplex(std::vector<double> &x, const std::vector<char> &active) {
    double sum = 0.0;
    for (size_t g = 0; g < x.size(); ++g) {
        if (!active[g] || x[g] < 0.0) {
            x[g] = 0.0;
        }
        sum += x[g];
    }
    if (sum > 0.0) {
        for (size_t g = 0; g < x.size(); ++g) {
            x[g] /= sum;
        }
    }
}

void runEm(const Model &m, EmState &s, const std::vector<double> &invLen, double priorMass, int maxIter, double tol,
           int threads) {
    const size_t U = m.units(), G = m.groupCount();
    size_t activeGroups = 0;
    for (size_t g = 0; g < G; ++g) {
        activeGroups += s.active[g] ? 1 : 0;
    }
    const double priorTerm = activeGroups > 0 ? std::max(0.0, priorMass) * static_cast<double>(U) / static_cast<double>(activeGroups) : 0.0;

    // start: every fragment split over its candidates in proportion to phi
    std::vector<double> x0(G, 0.0);
    for (size_t u = 0; u < U; ++u) {
        double den = 0.0;
        for (size_t i = m.ugOffset[u]; i < m.ugOffset[u + 1]; ++i) {
            den += s.active[m.ugGroup[i]] ? m.ugPhi[i] : 0.0;
        }
        for (size_t i = m.ugOffset[u]; den > 0.0 && i < m.ugOffset[u + 1]; ++i) {
            if (s.active[m.ugGroup[i]]) {
                x0[m.ugGroup[i]] += m.ugPhi[i] / den;
            }
        }
    }
    projectSimplex(x0, s.active);

    double stepMax = STEP_LENGTH_MIN;
    int iter = 0;
    bool converged = false;
    for (; iter < maxIter; ++iter) {
        const std::vector<double> x1 = emMap(m, s, x0, invLen, priorTerm, threads);
        const std::vector<double> x2 = emMap(m, s, x1, invLen, priorTerm, threads);
        double nr = 0.0, nv = 0.0;
        for (size_t g = 0; g < G; ++g) {
            const double r = x1[g] - x0[g];
            const double v = x2[g] - 2.0 * x1[g] + x0[g];
            nr += r * r;
            nv += v * v;
        }
        double a = nv > 0.0 ? -std::sqrt(nr / nv) : -STEP_LENGTH_MIN;
        if (!std::isfinite(a)) {
            a = -STEP_LENGTH_MIN;
        }
        a = std::max(-stepMax, std::min(-STEP_LENGTH_MIN, a));
        const bool capBinding = a <= -stepMax + 1e-12;
        std::vector<double> xn(G);
        for (size_t g = 0; g < G; ++g) {
            const double r = x1[g] - x0[g];
            const double v = x2[g] - 2.0 * x1[g] + x0[g];
            xn[g] = x0[g] - 2.0 * a * r + a * a * v;
        }
        projectSimplex(xn, s.active);
        const double llEm = logLikelihood(m, s, x2, invLen, threads);
        const double llAcc = logLikelihood(m, s, xn, invLen, threads);
        if (std::isfinite(llAcc) && llAcc >= llEm - 1e-12) {
            if (capBinding) {
                stepMax = std::min(STEP_LENGTH_CAP, stepMax * STEP_LENGTH_GROWTH);
            }
        } else {
            xn = x2;
            stepMax = std::max(STEP_LENGTH_MIN, stepMax / STEP_LENGTH_GROWTH);
        }
        double change = 0.0;
        for (size_t g = 0; g < G; ++g) {
            change = std::max(change, std::fabs(xn[g] - x0[g]));
        }
        x0.swap(xn);
        if (change < tol && iter + 1 >= MIN_ITERATIONS) {
            converged = true;
            ++iter;
            break;
        }
    }
    s.theta = x0;
    // posteriors at the final theta (the M-step output is not needed)
    emMap(m, s, s.theta, invLen, priorTerm, threads);
    Debug(Debug::INFO) << "EM over " << U << " fragments and " << activeGroups << " groups: "
                       << (converged ? "converged after " : "stopped at the iteration limit, ") << iter << " SQUAREM iterations.\n";
}

// gamma of the (fragment of read q, group of hit h) candidate, 0 if the pair rule removed it
double fragmentGroupPosterior(const Model &m, const EmState &s, unsigned int q, size_t h) {
    const unsigned int u = m.unitOfQuery[q];
    for (size_t i = m.ugOffset[u]; i < m.ugOffset[u + 1]; ++i) {
        if (m.ugGroup[i] == m.hitGroup[h]) {
            return s.gamma[i];
        }
    }
    return 0.0;
}

} // namespace

int emreclassify(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    Timer timer;

    // ---- genomes (--ref-group 1): groups and their lengths over the WHOLE target DB (hit or not)
    std::vector<unsigned int> genomeOfKey;
    std::vector<double> genomeLength;
    if (par.emRefGroup) {
        DBReader<unsigned int> tdb(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX);
        tdb.open(DBReader<unsigned int>::NOSORT);
        unsigned int maxKey = 0;
        for (size_t i = 0; i < tdb.getSize(); ++i) {
            maxKey = std::max(maxKey, tdb.getDbKey(i));
        }
        EMTargetGroups::Groups groups = EMTargetGroups::load(par.db2, par.emTargetGroups, maxKey);
        genomeLength.assign(groups.names.size(), 0.0);
        for (size_t i = 0; i < tdb.getSize(); ++i) {
            const unsigned int g = groups.groupOfKey[tdb.getDbKey(i)];
            if (g != NONE) {
                genomeLength[g] += static_cast<double>(std::max<size_t>(1, tdb.getSeqLen(i)));
            }
        }
        Debug(Debug::INFO) << "Genomes from " << groups.how << ": " << groups.mapped << " of " << tdb.getSize()
                           << " targets in " << groups.names.size() << " genomes, the others each on their own.\n";
        tdb.close();
        genomeOfKey.swap(groups.groupOfKey);
    }

    // ---- hits, fragments, compatibilities
    Model m;
    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    loadHits(reader, m, par.threads);
    assignGroups(m, genomeOfKey, genomeLength, par.threads);
    std::vector<unsigned int>().swap(genomeOfKey);
    Debug(Debug::INFO) << "Loaded " << m.queryKey.size() << " reads with " << m.queryOffset.back() << " hits.\n";
    buildFragments(m, par.db1, par.emMatePair != 0);
    buildCompatibilities(m, par.reclassifyLambda, par.threads);

    // ---- EM
    const size_t Q = m.queryKey.size(), G = m.groupCount();
    std::vector<double> invLen(G);
    for (size_t g = 0; g < G; ++g) {
        invLen[g] = 1.0 / std::max(1.0, m.groupLength[g]);
    }
    EmState s;
    s.gamma.assign(m.ugOffset.back(), 0.0);
    s.active.assign(G, 0);
    for (size_t i = 0; i < m.ugGroup.size(); ++i) {
        s.active[m.ugGroup[i]] = 1;               // genomes without a candidate stay out of the EM
    }
    runEm(m, s, invLen, par.reclassifyPriorMass, par.reclassifyMaxIterations, par.reclassifyTolerance, par.threads);

    // ---- read-level posteriors and output
    std::vector<double> hitPost(m.queryOffset.back(), 0.0);
#pragma omp parallel for schedule(dynamic, 256) num_threads(usableThreads(par.threads))
    for (size_t q = 0; q < Q; ++q) {
        double sum = 0.0;
        for (size_t h = m.queryOffset[q]; h < m.queryOffset[q + 1]; ++h) {
            hitPost[h] = fragmentGroupPosterior(m, s, static_cast<unsigned int>(q), h) * m.hitPhi[h] / m.hitPhiQG[h];
            sum += hitPost[h];
        }
        for (size_t h = m.queryOffset[q]; h < m.queryOffset[q + 1]; ++h) {
            hitPost[h] = sum > 0.0 ? hitPost[h] / sum : 0.0;
        }
    }

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();
    size_t assigned = 0;
#pragma omp parallel num_threads(usableThreads(par.threads)) reduction(+:assigned)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024 + 32768 * 4];
        std::vector<size_t> order;
#pragma omp for schedule(dynamic, 16)
        for (size_t q = 0; q < Q; ++q) {
            order.clear();
            for (size_t h = m.queryOffset[q]; h < m.queryOffset[q + 1]; ++h) {
                if (hitPost[h] > 0.0) {
                    order.push_back(h);
                }
            }
            SORT_SERIAL(order.begin(), order.end(), [&](size_t a, size_t b) {
                if (hitPost[a] != hitPost[b]) {
                    return hitPost[a] > hitPost[b];
                }
                if (m.hitScore[a] != m.hitScore[b]) {
                    return m.hitScore[a] > m.hitScore[b];
                }
                return Matcher::compareHits(m.hitResult[a], m.hitResult[b]);
            });
            assigned += order.empty() ? 0 : 1;
            writer.writeStart(thread_idx);
            for (size_t j = 0; j < order.size(); ++j) {
                Matcher::result_t r = m.hitResult[order[j]];
                r.seqId = static_cast<float>(hitPost[order[j]]);
                const size_t len = Matcher::resultToBuffer(buffer, r, m.hasBacktrace, false, m.hasOrfPosition, true);
                writer.writeAdd(buffer, len, thread_idx);
            }
            writer.writeEnd(m.queryKey[q], thread_idx);
        }
    }
    writer.close();
    reader.close();
    Debug(Debug::INFO) << assigned << " of " << Q << " reads assigned. Time: " << timer.lap() << "\n";
    return EXIT_SUCCESS;
}
