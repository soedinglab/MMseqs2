// Resource monitor for the one-shot workflow. Samples /proc only, so it runs unchanged on any
// Linux machine, and follows the process tree below one pid: the workflow shell and every
// mmseqs pass it starts. One TSV line per interval, a per-phase summary when the tree is gone.
#include "Parameters.h"
#include "Debug.h"
#include "Util.h"

#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <dirent.h>
#include <map>
#include <set>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {

volatile sig_atomic_t stopRequested = 0;
void onStop(int) { stopRequested = 1; }

struct ProcessSample {
    pid_t parent;
    std::string command;
    unsigned long long cpuTicks;
    unsigned long long rssBytes;
    unsigned long long readChars;
    unsigned long long writeChars;
    unsigned long long readBytes;
    unsigned long long writeBytes;
};

struct NodeSample {
    unsigned long long cpuBusy;
    unsigned long long cpuTotal;
    unsigned long long memUsed;
    unsigned long long memTotal;
    unsigned long long netRx;
    unsigned long long netTx;
};

struct PhaseTotals {
    double seconds;
    double cpuSeconds;
    unsigned long long peakRss;
    unsigned long long readChars;
    unsigned long long writeChars;
    PhaseTotals() : seconds(0), cpuSeconds(0), peakRss(0), readChars(0), writeChars(0) {}
};

bool readWhole(const std::string &path, std::string &into) {
    FILE *in = fopen(path.c_str(), "r");
    if (in == NULL) {
        return false;
    }
    char buffer[4096];
    into.clear();
    size_t got;
    while ((got = fread(buffer, 1, sizeof(buffer), in)) > 0) {
        into.append(buffer, got);
    }
    fclose(in);
    return true;
}

bool readNode(NodeSample &node) {
    std::string text;
    if (readWhole("/proc/stat", text) == false) {
        return false;
    }
    unsigned long long v[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    if (sscanf(text.c_str(), "cpu %llu %llu %llu %llu %llu %llu %llu %llu", &v[0], &v[1], &v[2], &v[3],
               &v[4], &v[5], &v[6], &v[7]) < 4) {
        return false;
    }
    node.cpuTotal = 0;
    for (int i = 0; i < 8; i++) {
        node.cpuTotal += v[i];
    }
    node.cpuBusy = node.cpuTotal - v[3] - v[4];

    node.memTotal = 0;
    unsigned long long available = 0;
    if (readWhole("/proc/meminfo", text)) {
        const char *at = strstr(text.c_str(), "MemTotal:");
        if (at != NULL) {
            sscanf(at, "MemTotal: %llu", &node.memTotal);
        }
        at = strstr(text.c_str(), "MemAvailable:");
        if (at != NULL) {
            sscanf(at, "MemAvailable: %llu", &available);
        }
    }
    node.memTotal *= 1024;
    node.memUsed = node.memTotal - available * 1024;

    node.netRx = 0;
    node.netTx = 0;
    if (readWhole("/proc/net/dev", text)) {
        size_t line = 0;
        while (line < text.size()) {
            size_t end = text.find('\n', line);
            if (end == std::string::npos) {
                end = text.size();
            }
            const std::string row = text.substr(line, end - line);
            line = end + 1;
            const size_t colon = row.find(':');
            if (colon == std::string::npos) {
                continue;
            }
            std::string name = row.substr(0, colon);
            name.erase(0, name.find_first_not_of(" \t"));
            if (name == "lo") {
                continue;
            }
            unsigned long long rx = 0, tx = 0, skip = 0;
            if (sscanf(row.c_str() + colon + 1, "%llu %llu %llu %llu %llu %llu %llu %llu %llu", &rx, &skip, &skip,
                       &skip, &skip, &skip, &skip, &skip, &tx) == 9) {
                node.netRx += rx;
                node.netTx += tx;
            }
        }
    }
    return true;
}

bool readProcess(pid_t pid, ProcessSample &sample) {
    std::string text;
    if (readWhole("/proc/" + SSTR(pid) + "/stat", text) == false) {
        return false;
    }
    // the command sits in parentheses and may hold spaces, so cut at the last one
    const size_t open = text.find('(');
    const size_t close = text.rfind(')');
    if (open == std::string::npos || close == std::string::npos || close < open) {
        return false;
    }
    sample.command = text.substr(open + 1, close - open - 1);
    const char *rest = text.c_str() + close + 2;
    char state = 0;
    long ppid = 0;
    unsigned long long utime = 0, stime = 0;
    long rssPages = 0;
    // fields 3.. of the stat line: state ppid pgrp session tty tpgid flags minflt cminflt majflt
    // cmajflt utime stime cutime cstime priority nice threads itreal start vsize rss
    if (sscanf(rest, "%c %ld %*d %*d %*d %*d %*u %*u %*u %*u %*u %llu %llu %*d %*d %*d %*d %*d %*d %*u %*u %ld",
               &state, &ppid, &utime, &stime, &rssPages) != 5) {
        return false;
    }
    sample.parent = (pid_t) ppid;
    sample.cpuTicks = utime + stime;
    sample.rssBytes = (unsigned long long) rssPages * (unsigned long long) sysconf(_SC_PAGESIZE);
    sample.readChars = sample.writeChars = sample.readBytes = sample.writeBytes = 0;
    if (readWhole("/proc/" + SSTR(pid) + "/io", text)) {
        const char *at;
        if ((at = strstr(text.c_str(), "rchar:")) != NULL) sscanf(at, "rchar: %llu", &sample.readChars);
        if ((at = strstr(text.c_str(), "wchar:")) != NULL) sscanf(at, "wchar: %llu", &sample.writeChars);
        if ((at = strstr(text.c_str(), "read_bytes:")) != NULL) sscanf(at, "read_bytes: %llu", &sample.readBytes);
        if ((at = strstr(text.c_str(), "write_bytes:")) != NULL) sscanf(at, "write_bytes: %llu", &sample.writeBytes);
    }
    return true;
}

// the mmseqs pass a process runs, from its argv[1]; other programs keep their name
std::string phaseName(pid_t pid, const std::string &command) {
    if (command != "mmseqs") {
        return command;
    }
    std::string text;
    if (readWhole("/proc/" + SSTR(pid) + "/cmdline", text) == false) {
        return command;
    }
    size_t at = text.find('\0');
    if (at == std::string::npos || at + 1 >= text.size()) {
        return "";
    }
    const size_t end = text.find('\0', at + 1);
    return text.substr(at + 1, end == std::string::npos ? std::string::npos : end - at - 1);
}

// every process below root, found by one scan of /proc
void readTree(pid_t root, std::map<pid_t, ProcessSample> &tree, std::set<std::string> &phases) {
    tree.clear();
    phases.clear();
    const pid_t self = getpid();
    std::string rootName;
    std::map<pid_t, ProcessSample> all;
    DIR *dir = opendir("/proc");
    if (dir == NULL) {
        return;
    }
    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[0] < '0' || entry->d_name[0] > '9') {
            continue;
        }
        const pid_t pid = (pid_t) atol(entry->d_name);
        ProcessSample sample;
        if (readProcess(pid, sample)) {
            all[pid] = sample;
        }
    }
    closedir(dir);
    std::map<pid_t, std::vector<pid_t> > children;
    for (std::map<pid_t, ProcessSample>::const_iterator it = all.begin(); it != all.end(); ++it) {
        children[it->second.parent].push_back(it->first);
    }
    if (all.count(root) != 0) {
        rootName = all[root].command;
    }
    std::vector<pid_t> queue(1, root);
    while (queue.empty() == false) {
        const pid_t pid = queue.back();
        queue.pop_back();
        const std::map<pid_t, ProcessSample>::const_iterator me = all.find(pid);
        if (me == all.end()) {
            continue;
        }
        if (pid == self) {
            continue;
        }
        tree[pid] = me->second;
        const std::string &name = me->second.command;
        if (pid != root && name != rootName && name != "sh" && name != "bash" && name != "sleep") {
            const std::string phase = phaseName(pid, name);
            if (phase.empty() == false) {
                phases.insert(phase);
            }
        }
        const std::map<pid_t, std::vector<pid_t> >::const_iterator kids = children.find(pid);
        if (kids != children.end()) {
            queue.insert(queue.end(), kids->second.begin(), kids->second.end());
        }
    }
}

bool alive(pid_t pid) {
    return kill(pid, 0) == 0 || errno == EPERM;
}

std::string isoTime(time_t at) {
    char buffer[32];
    struct tm parts;
    localtime_r(&at, &parts);
    strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S", &parts);
    return buffer;
}

double gb(unsigned long long bytes) { return (double) bytes / (1024.0 * 1024.0 * 1024.0); }
double mbPerSecond(unsigned long long bytes, double seconds) {
    return seconds > 0 ? (double) bytes / (1024.0 * 1024.0) / seconds : 0;
}

}

int lin8monitor(int argc, const char **argv, const Command &command) {
    signal(SIGTERM, onStop);
    signal(SIGINT, onStop);
    signal(SIGHUP, onStop);
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
#ifndef __linux__
    Debug(Debug::ERROR) << "lin8-monitor reads /proc and runs on Linux only\n";
    return EXIT_FAILURE;
#endif
    if (par.lin8MonitorPid <= 0) {
        Debug(Debug::ERROR) << "--monitor-pid must name the process whose tree to watch\n";
        return EXIT_FAILURE;
    }
    const pid_t root = (pid_t) par.lin8MonitorPid;
    const int interval = std::max(1, par.lin8MonitorInterval);
    const double ticksPerSecond = (double) sysconf(_SC_CLK_TCK);
    const long cpus = sysconf(_SC_NPROCESSORS_ONLN);

    struct stat sb;
    const bool fresh = stat(par.db1.c_str(), &sb) != 0 || sb.st_size == 0;
    FILE *out = fopen(par.db1.c_str(), "a");
    if (out == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << par.db1 << " for appending\n";
        return EXIT_FAILURE;
    }
    setvbuf(out, NULL, _IOLBF, 1 << 16);

    char host[256];
    memset(host, 0, sizeof(host));
    gethostname(host, sizeof(host) - 1);
    const time_t started = time(NULL);
    NodeSample node;
    readNode(node);
    fprintf(out, "# lin8-monitor host=%s cpus=%ld mem_total_gb=%.1f pid=%d interval_s=%d %s=%s\n", host, cpus,
            gb(node.memTotal), (int) root, interval, fresh ? "started" : "resumed", isoTime(started).c_str());
    if (fresh) {
        fprintf(out, "time\telapsed_s\tphase\tnode_cpu_cores\ttree_cpu_cores\tnode_mem_used_gb\ttree_rss_gb\t"
                     "tree_rss_peak_gb\ttree_read_mb_s\ttree_write_mb_s\tdisk_read_mb_s\tdisk_write_mb_s\t"
                     "net_rx_mb_s\tnet_tx_mb_s\n");
    }

    NodeSample previousNode = node;
    std::map<pid_t, ProcessSample> previousTree;
    std::set<std::string> phases;
    readTree(root, previousTree, phases);
    time_t previousTime = started;
    unsigned long long peakRss = 0;
    double totalCpuSeconds = 0;
    unsigned long long totalRead = 0, totalWrite = 0, totalDiskRead = 0, totalDiskWrite = 0;
    unsigned long long totalNetRx = 0, totalNetTx = 0;
    std::map<std::string, PhaseTotals> perPhase;
    std::vector<std::string> phaseOrder;

    while (stopRequested == 0) {
        for (int slept = 0; slept < interval && stopRequested == 0 && alive(root); slept++) {
            sleep(1);
        }
        const bool rootAlive = alive(root);
        const time_t now = time(NULL);
        const double seconds = (double) (now - previousTime);
        std::map<pid_t, ProcessSample> tree;
        readTree(root, tree, phases);
        readNode(node);

        unsigned long long rss = 0, cpuTicks = 0, readChars = 0, writeChars = 0, readBytes = 0, writeBytes = 0;
        for (std::map<pid_t, ProcessSample>::const_iterator it = tree.begin(); it != tree.end(); ++it) {
            const ProcessSample &s = it->second;
            rss += s.rssBytes;
            const std::map<pid_t, ProcessSample>::const_iterator before = previousTree.find(it->first);
            // a process seen for the first time started after the last sample, so all of it counts
            const ProcessSample zero = {0, "", 0, 0, 0, 0, 0, 0};
            const ProcessSample &p = (before == previousTree.end()) ? zero : before->second;
            cpuTicks += s.cpuTicks - std::min(s.cpuTicks, p.cpuTicks);
            readChars += s.readChars - std::min(s.readChars, p.readChars);
            writeChars += s.writeChars - std::min(s.writeChars, p.writeChars);
            readBytes += s.readBytes - std::min(s.readBytes, p.readBytes);
            writeBytes += s.writeBytes - std::min(s.writeBytes, p.writeBytes);
        }
        peakRss = std::max(peakRss, rss);
        const double cpuSeconds = (double) cpuTicks / ticksPerSecond;
        totalCpuSeconds += cpuSeconds;
        totalRead += readChars;
        totalWrite += writeChars;
        totalDiskRead += readBytes;
        totalDiskWrite += writeBytes;
        const unsigned long long netRx = node.netRx - std::min(node.netRx, previousNode.netRx);
        const unsigned long long netTx = node.netTx - std::min(node.netTx, previousNode.netTx);
        totalNetRx += netRx;
        totalNetTx += netTx;
        const unsigned long long busy = node.cpuBusy - std::min(node.cpuBusy, previousNode.cpuBusy);
        const unsigned long long total = node.cpuTotal - std::min(node.cpuTotal, previousNode.cpuTotal);
        const double nodeCores = total > 0 ? (double) busy / (double) total * (double) cpus : 0;

        std::string phase;
        for (std::set<std::string>::const_iterator it = phases.begin(); it != phases.end(); ++it) {
            phase += (phase.empty() ? "" : "+") + *it;
        }
        if (phase.empty()) {
            phase = rootAlive ? "-" : "exit";
        }
        if (perPhase.count(phase) == 0) {
            phaseOrder.push_back(phase);
        }
        PhaseTotals &totals = perPhase[phase];
        totals.seconds += seconds;
        totals.cpuSeconds += cpuSeconds;
        totals.peakRss = std::max(totals.peakRss, rss);
        totals.readChars += readChars;
        totals.writeChars += writeChars;

        fprintf(out, "%s\t%ld\t%s\t%.2f\t%.2f\t%.2f\t%.3f\t%.3f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",
                isoTime(now).c_str(), (long) (now - started), phase.c_str(), nodeCores,
                seconds > 0 ? cpuSeconds / seconds : 0, gb(node.memUsed), gb(rss), gb(peakRss),
                mbPerSecond(readChars, seconds), mbPerSecond(writeChars, seconds),
                mbPerSecond(readBytes, seconds), mbPerSecond(writeBytes, seconds),
                mbPerSecond(netRx, seconds), mbPerSecond(netTx, seconds));
        previousTree.swap(tree);
        previousNode = node;
        previousTime = now;
        if (rootAlive == false) {
            break;
        }
    }

    const time_t ended = time(NULL);
    fprintf(out, "# summary elapsed_s=%ld tree_cpu_s=%.0f tree_cpu_cores_avg=%.2f tree_rss_peak_gb=%.3f "
                 "tree_read_gb=%.2f tree_write_gb=%.2f disk_read_gb=%.2f disk_write_gb=%.2f "
                 "net_rx_gb=%.2f net_tx_gb=%.2f ended=%s\n",
            (long) (ended - started), totalCpuSeconds,
            ended > started ? totalCpuSeconds / (double) (ended - started) : 0, gb(peakRss), gb(totalRead),
            gb(totalWrite), gb(totalDiskRead), gb(totalDiskWrite), gb(totalNetRx), gb(totalNetTx),
            isoTime(ended).c_str());
    for (size_t i = 0; i < phaseOrder.size(); i++) {
        const PhaseTotals &t = perPhase[phaseOrder[i]];
        fprintf(out, "# phase %s seconds=%.0f cpu_s=%.0f cpu_cores_avg=%.2f rss_peak_gb=%.3f read_gb=%.2f write_gb=%.2f\n",
                phaseOrder[i].c_str(), t.seconds, t.cpuSeconds, t.seconds > 0 ? t.cpuSeconds / t.seconds : 0,
                gb(t.peakRss), gb(t.readChars), gb(t.writeChars));
    }
    fclose(out);
    return EXIT_SUCCESS;
}
