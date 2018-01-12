#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"

#include <algorithm>
#include <climits>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <spawn.h>
#include <poll.h>

#ifdef OPENMP
#include <omp.h>
#endif

// Analogous to gnulib implementation
// https://www.gnu.org/software/gnulib/
// Licensed under GPLv3
pid_t create_pipe(const char *prog_path, char **prog_argv, char **environ, int fd[2]) {
    int ifd[2];
    int ofd[2];

    int err;
    if ((err = pipe(ifd)) < 0) {
        perror("pipe ifd");
        errno = err;
        return -1;
    }
    if ((err = pipe(ofd)) < 0) {
        perror("pipe ofd");
        errno = err;
        return -1;
    }

    int actions_allocated = 0;
    int attrs_allocated = 0;

    posix_spawn_file_actions_t actions;
    posix_spawnattr_t attrs;
    pid_t child;

    if ((err = posix_spawn_file_actions_init(&actions)) != 0
        || (actions_allocated = 1,
            (err = posix_spawn_file_actions_adddup2 (&actions, ofd[0], STDIN_FILENO)) != 0
            || (err = posix_spawn_file_actions_adddup2 (&actions, ifd[1], STDOUT_FILENO)) != 0
            || (err = posix_spawn_file_actions_addclose(&actions, ofd[0])) != 0
            || (err = posix_spawn_file_actions_addclose(&actions, ifd[1])) != 0
            || (err = posix_spawn_file_actions_addclose(&actions, ofd[1])) != 0
            || (err = posix_spawn_file_actions_addclose(&actions, ifd[0])) != 0
            #ifdef POSIX_SPAWN_USEVFORK
            || ((err = posix_spawnattr_init(&attrs)) != 0
                || (attrs_allocated = 1,
                   (err = posix_spawnattr_setflags(&attrs, POSIX_SPAWN_USEVFORK)) != 0))
            #endif
            || (err = posix_spawnp(&child, prog_path, &actions, attrs_allocated ? &attrs : NULL, prog_argv, environ)) != 0))
    {
        perror("fail");
        errno = err;

        if (actions_allocated) {
            posix_spawn_file_actions_destroy(&actions);
        }
        if (attrs_allocated) {
            posix_spawnattr_destroy(&attrs);
        }

        close(ifd[0]);
        close(ifd[1]);
        close(ofd[0]);
        close(ofd[1]);

        return -1;
    }

    posix_spawn_file_actions_destroy(&actions);
    if (attrs_allocated) {
        posix_spawnattr_destroy(&attrs);
    }

    if ((err = close(ofd[0])) == -1 || (err = close(ifd[1])) == -1) {
        perror("close");
        errno = err;
        return -1;
    }

    fd[0] = ifd[0];
    fd[1] = ofd[1];
    return child;
}

int apply_by_entry(char* data, size_t size, unsigned int key, DBWriter& writer,
                   const char* program_name, char ** program_argv, char **environ, unsigned int thread_idx) {
    // only works with the environ we construct ourselves
    // local_environment() leaves the first element free to use for ourselves
    snprintf(environ[0], 64, "MMSEQS_ENTRY_NAME=%d", key);

    bool write_closed = false;
    int fd[2];
    pid_t child_pid;
    if ((child_pid = create_pipe(program_name, program_argv, environ, fd)) == -1) {
        perror("create_pipe");
        return -1;
    }

    // Analogous to gnulib implementation
    size_t written = 0;
    int error = 0;
    int fcntl_flags;
    if ( (fcntl_flags = fcntl(fd[1], F_GETFL, 0)) < 0
         || fcntl(fd[1], F_SETFL, fcntl_flags | O_NONBLOCK) == -1
         || (fcntl_flags = fcntl(fd[0], F_GETFL, 0)) < 0
         || fcntl(fd[0], F_SETFL, fcntl_flags | O_NONBLOCK) == -1)
    {
        goto end;
    }

    char buffer[PIPE_BUF];
    writer.writeStart(thread_idx);

    struct pollfd plist[2];
    for (;;) {
        size_t rest = size - written;
        size_t batch_size = PIPE_BUF;
        if (rest < PIPE_BUF) {
            batch_size = rest;
        }

        plist[0].fd = write_closed == false ? fd[1] : fd[1] * -1;
        plist[0].events = POLLOUT;
        plist[0].revents = 0;

        plist[1].fd = fd[0];
        plist[1].events = POLLIN;
        plist[1].revents = 0;

        if (poll(plist, 2, -1) == -1) {
            if (errno == EAGAIN) {
                perror("again");
                continue;
            }
            perror("poll");
            error = errno;
            break;
        }

        if (plist[0].revents & POLLOUT) {
            ssize_t write_size = batch_size;
            if (size - written > 0) {
                for(;;) {
                    ssize_t w = write(fd[1], data + written, write_size);
                    if (w < 0) {
                        if (errno != EAGAIN) {
                            perror("write stdin1");
                            error = errno;
                            goto end;
                        } else {
                            write_size = write_size / 2;
                            if (write_size == 0) {
                                break;
                            }
                        }
                    } else {
                        written += w;
                        break;
                    }
                }
            } else {
                if (close(fd[1]) == -1) {
                    perror("close error");
                    error = errno;
                    break;
                }
                write_closed = true;
            }
        } else if (plist[1].revents & POLLIN) {
            ssize_t bytes_read = read(plist[1].fd, &buffer, sizeof(buffer));
            if (bytes_read > 0) {
                writer.writeAdd(buffer, bytes_read, thread_idx);
            } else if (bytes_read < 0) {
                if (errno != EAGAIN) {
                    perror("read stdout0");
                    error = errno;
                    break;
                }
            } else if (bytes_read == 0 && write_closed == true) {
                break;
            }
        } else {
            // nothing left to read or write
            break;
        }
    }

    end:

    writer.writeEnd(key, thread_idx, true);

    if (write_closed == true) {
        close(fd[1]);
    }

    if (close(fd[0]) == -1) {
        perror("close stdout");
        error = errno;
    }

    int status = 0;
    while (waitpid(child_pid, &status, 0) == -1) {
        if (errno == EINTR) {
            continue;
        }
        perror("waitpid");
        error = errno;
    }

    errno = error;
    return WEXITSTATUS(status);
}

void ignore_signal(int signal) {
    struct sigaction handler;
    handler.sa_handler = SIG_IGN;
    sigemptyset(&handler.sa_mask);
    handler.sa_flags = 0;
    sigaction(signal, &handler, NULL);
}

extern char **environ;

int prefix(const char *pre, const char *str) {
    return strncmp(pre, str, strlen(pre)) == 0;
}

char** local_environment() {
    size_t env_size = 0;
    for (size_t i = 0; environ[i] != NULL; ++i) {
        env_size++;
    }

    char** local_environ = (char**)malloc(sizeof(char*) * (env_size + 2));
    size_t j = 0;
    local_environ[j++] = (char*)malloc(sizeof(char) * 64);
    for (size_t i = 0; environ[i] != NULL; ++i) {
        if (prefix("MMSEQS_ENTRY_NAME=", environ[i]) == 0) {
            local_environ[j++] = environ[i];
        }
    }
    local_environ[j] = NULL;

    return local_environ;
}

void free_local_environment(char** local_environ) {
    free(local_environ[0]);
    free(local_environ);
}

int apply(int argc, const char **argv, const Command& command) {
    ignore_signal(SIGPIPE);

    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2, true, Parameters::PARSE_REST);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());

    // to make processing time as equal as possible over the different MPI processes and threads,
    // we shuffle the input database and later sort by entry length per thread
    // with the assumption that longer input will take longer to process and
    // that, if the most time consuming entries are processed first, no long stragglers in the end will
    // force the other processes / threads to stall
    reader.open(DBReader<unsigned int>::SHUFFLE);

#ifdef HAVE_MPI
    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(reader.getAminoAcidDBSize(), reader.getSeqLens(), reader.getSize(),
                                     MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db2, par.db2Index, MMseqsMPI::rank);
    const char* resultData = tmpOutput.first.c_str();
    const char* resultIndex = tmpOutput.second.c_str();
#else
    size_t dbFrom = 0;
    size_t dbSize = reader.getSize();
    const char* resultData = par.db2.c_str();
    const char* resultIndex = par.db2Index.c_str();
#endif

    DBWriter writer(resultData, resultIndex, par.threads);
    writer.open();

    unsigned int *sizes = reader.getSeqLens();
    Debug(Debug::INFO) << "Start applying.\n";
#pragma omp parallel
    {
        char **local_environ = local_environment();

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        std::pair<unsigned int, unsigned int> *lengthSorted = new std::pair<unsigned int, unsigned int>[dbSize];
        for (size_t i = dbFrom; i < (dbFrom + dbSize); i++) {
            lengthSorted[i - dbFrom] = std::make_pair(i, sizes[i]);
        }
        std::stable_sort(lengthSorted, lengthSorted + dbSize, DBReader<unsigned int>::comparePairBySeqLength());

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < dbSize; ++i) {
            Debug::printProgress(i);

            size_t index = lengthSorted[i].first;
            size_t size = lengthSorted[i].second - 1;

            char *data = reader.getData(index);
            if (data == NULL) {
                continue;
            }

            unsigned int key = reader.getDbKey(index);
            int status = apply_by_entry(data, size, key, writer, par.restArgv[0], const_cast<char**>(par.restArgv), local_environ, thread_idx);
            if (status == -1) {
                Debug(Debug::WARNING) << "Entry " << index << " system error " << errno << "!\n";
                continue;
            } else if (status > 0) {
                Debug(Debug::WARNING) << "Entry " << index << " exited with error code " << status << "!\n";
                continue;
            }
        }

        delete[] lengthSorted;

        free_local_environment(local_environ);
    }
    Debug(Debug::INFO) << "\n";
    writer.close();

    Debug(Debug::INFO) << "\nDone.\n";
    reader.close();

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    if(MMseqsMPI::rank == MMseqsMPI::MASTER) {
        // master reduces results
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (int proc = 0; proc < MMseqsMPI::numProc; ++proc) {
            splitFiles.emplace_back(Util::createTmpFileNames(par.db2, par.db2Index, proc));
        }
        DBWriter::mergeResults(par.db2, par.db2Index, splitFiles);
    }
#endif

    return EXIT_SUCCESS;
}

