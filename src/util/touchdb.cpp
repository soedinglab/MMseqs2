#include "Parameters.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "MemoryMapped.h"

#include <sys/mman.h>

int touchdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const char* mlockError = "Consider granting the `CAP_IPC_LOCK` capability to this binary (e.g., `setcap cap_ipc_lock=+ep <program>`, or run as root) or raising the locked-memory limit so `mlock(2)` can succeed. Check your current limit with `ulimit -l` and increase it via `ulimit -l <KB>` for the shell, `LimitMEMLOCK=infinity` in your systemd unit, or persistent user limits in `/etc/security/limits.conf` (e.g., `youruser hard memlock unlimited`)\n";

    std::string db = par.db1;
    std::string indexDB = PrefilteringIndexReader::searchForIndex(db);
    if (indexDB.empty() == false) {
        db = indexDB;
        if (par.idList != "") {
            std::string idx = db + ".index";
            DBReader<unsigned int> reader(db.c_str(), idx.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
            reader.open(DBReader<unsigned int>::NOSORT);

            std::vector<std::string> ids = Util::split(par.idList, ",");
            for (size_t i = 0; i < ids.size(); ++i) {
                size_t id = reader.getId(Util::fast_atoi<unsigned int>(ids[i].c_str()));
                if (id == UINT_MAX) {
                    Debug(Debug::WARNING) << "Key " << ids[i] << " not found in database\n";
                    continue;
                }
                size_t currDataOffset = reader.getOffset(id);
                size_t nextDataOffset = reader.findNextOffsetid(id);
                size_t dataSize = nextDataOffset - currDataOffset;
                char* data = reader.getDataUncompressed(id);
                Util::touchMemory(data, dataSize);
                if (par.touchLock) {
                    int res = mlock(data, dataSize);
                    if (res != 0) {
                        Debug(Debug::ERROR) << "Could not lock memory " << strerror(errno) << "\n";
                        Debug(Debug::ERROR) << mlockError;
                        reader.close();
                        return EXIT_FAILURE;
                    }
                }
            }

            if (par.touchLock) {
                Debug(Debug::INFO) << "Touched and locked " << ids.size() << " entries. Process will not exit until killed.\n";
                pause();
            } else {
                Debug(Debug::INFO) << "Touched " << ids.size() << " entries\n";
            }

            reader.close();
            return EXIT_SUCCESS;
        }
    }

    MemoryMapped map(db, MemoryMapped::WholeFile, MemoryMapped::CacheHint::SequentialScan);
    Util::touchMemory(reinterpret_cast<const char*>(map.getData()), map.mappedSize());
    if (par.touchLock) {
        int res = mlock(map.getData(), map.mappedSize());
        if (res != 0) {
            Debug(Debug::ERROR) << "Could not lock memory " << strerror(errno) << "\n";
            Debug(Debug::ERROR) << mlockError;
            return EXIT_FAILURE;
        }
        Debug(Debug::INFO) << "Touched and locked database. Process will not exit until killed.\n";
        pause();
    }

    return EXIT_SUCCESS;
}
