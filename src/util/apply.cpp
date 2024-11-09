#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"

#if defined(__CYGWIN__) || defined(__EMSCRIPTEN__)
int apply(int, const char **, const Command&) {
    Debug(Debug::ERROR) << "\"apply\" is not supported on this platform\n";
    EXIT(EXIT_FAILURE);
}
#else
#include <climits>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <spawn.h>
#include <poll.h>
#include <sys/wait.h>
#include <sys/mman.h>

#ifdef OPENMP
#include <omp.h>
#endif

pid_t create_pipe(
    const char* prog_path,
    char** prog_argv,
    char** local_environ,
    int fd[2]) {
    if (prog_path == NULL) {
        return -1;
    }
    if (prog_argv == NULL) {
        return -1;
    }

    int pipe1Ids[2];
    int res = pipe(pipe1Ids);
    if (res != 0) {
        perror("pipe failed");
        return -1;
    }

    if (fcntl(pipe1Ids[0], F_SETFD, FD_CLOEXEC) == -1 || fcntl(pipe1Ids[1], F_SETFD, FD_CLOEXEC) == -1) {
        perror("fcntl failed");
        close(pipe1Ids[0]);
        close(pipe1Ids[1]);
        return -1;
    }

    int pipe2Ids[2];
    res = pipe(pipe2Ids);
    if (res != 0) {
        perror("pipe failed");
        close(pipe1Ids[0]);
        close(pipe1Ids[1]);
        return -1;
    }

    if (fcntl(pipe2Ids[0], F_SETFD, FD_CLOEXEC) == -1 || fcntl(pipe2Ids[1], F_SETFD, FD_CLOEXEC) == -1) {
        perror("fcntl failed");
        close(pipe1Ids[0]);
        close(pipe1Ids[1]);
        close(pipe2Ids[0]);
        close(pipe2Ids[1]);
        return -1;
    }

    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    if (posix_spawn_file_actions_adddup2(&actions, pipe2Ids[0], STDIN_FILENO) != 0
     || posix_spawn_file_actions_adddup2(&actions, pipe1Ids[1], STDOUT_FILENO) != 0) {
        perror("posix_spawn_file_actions failed");
        posix_spawn_file_actions_destroy(&actions);
        close(pipe1Ids[0]);
        close(pipe1Ids[1]);
        close(pipe2Ids[0]);
        close(pipe2Ids[1]);
        return -1;
     }

    int pid;
    res = posix_spawnp(&pid, prog_path, &actions, NULL, prog_argv, local_environ);
    if (res != 0) {
        perror("posix_spawn failed");
        close(pipe1Ids[0]);
        close(pipe1Ids[1]);
        close(pipe2Ids[0]);
        close(pipe2Ids[1]);
        return -1;
    }
    posix_spawn_file_actions_destroy(&actions);

    if (close(pipe2Ids[0]) != 0 || close(pipe1Ids[1]) != 0) {
        perror("close");
        return -1;
    }

    fd[0] = pipe1Ids[0];
    fd[1] = pipe2Ids[1];

    return pid;
}

int apply_by_entry(char* data, size_t size, unsigned int key, DBWriter& writer,
                   const char* program_name, char ** program_argv, char **environ, unsigned int proc_idx) {
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

    size_t written = 0;
    int error = 0;
    char buffer[PIPE_BUF];
    writer.writeStart(proc_idx);
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
                for (;;) {
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
                writer.writeAdd(buffer, bytes_read, proc_idx);
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

    writer.writeEnd(key, proc_idx, true);

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
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_REST, 0);

#ifdef OPENMP
    // forking does not play well with OpenMP threads
    omp_set_num_threads(1);
#endif

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::SORT_BY_LENGTH);

    Debug::Progress progress(reader.getSize());

#ifdef HAVE_MPI
    struct worker_s {
        int ready;
        int mpiRank;
        int mpiProc;
    };

    worker_s* shared_memory = (worker_s*)mmap(NULL, sizeof(worker_s), PROT_READ|PROT_WRITE, MAP_SHARED|MAP_ANONYMOUS, -1, 0);
    memset(shared_memory, 0, sizeof(worker_s));
#endif

    Debug(Debug::INFO) << "Start applying.\n";
    for (int thread = 0; thread < par.threads; ++thread) {
        switch (fork()) {
            default:
                break;
            case -1:
                Debug(Debug::ERROR) << "Could not fork worker process!\n";
                EXIT(EXIT_FAILURE);
            case 0: {
                int mpiRank = 0;
                int mpiProcs = 1;

                // get progress only from first thread on master
                if (thread != 0 || mpiRank != MMseqsMPI::MASTER) {
                    Debug::setDebugLevel(0);
                }

#ifdef HAVE_MPI
                while (shared_memory->ready == 0) {
                    usleep(10);
                }

                mpiRank = shared_memory->mpiRank;
                mpiProcs = shared_memory->mpiProc;
                std::pair<std::string, std::string> mpiDb = Util::createTmpFileNames(par.db2, par.db2Index, shared_memory->mpiRank);
                std::pair<std::string, std::string> outDb = Util::createTmpFileNames(mpiDb.first, mpiDb.second, thread);
#else
                std::pair<std::string, std::string> outDb = Util::createTmpFileNames(par.db2, par.db2Index, thread);
#endif

                DBWriter writer(outDb.first.c_str(), outDb.second.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
                writer.open();

                char **local_environ = local_environment();

                ignore_signal(SIGPIPE);
                for (size_t i = 0; i < reader.getSize(); ++i) {
                    progress.updateProgress();
                    if (static_cast<ssize_t>(i) % (mpiProcs * par.threads) != (thread * mpiProcs + mpiRank)) {
                        continue;
                    }

                    unsigned int key = reader.getDbKey(i);
                    char *data = reader.getData(i, thread);
                    if (*data == '\0') {
                        writer.writeData(NULL, 0, key, 0);
                        continue;
                    }

                    size_t size = reader.getEntryLen(i) - 1;
                    int status = apply_by_entry(data, size, key, writer, par.restArgv[0], const_cast<char**>(par.restArgv), local_environ, 0);
                    if (status == -1) {
                        Debug(Debug::WARNING) << "Entry " << key << " system error number " << errno << "!\n";
                        continue;
                    }
                    if (status > 0) {
                        Debug(Debug::WARNING) << "Entry " << key << " exited with error code " << status << "!\n";
                        continue;
                    }
                }

                writer.close(true);
                reader.close();
                free_local_environment(local_environ);
                _Exit(0);
            }
        }
    }

#ifdef HAVE_MPI
    shared_memory->mpiRank = MMseqsMPI::rank;
    shared_memory->mpiProc = MMseqsMPI::numProc;
    __sync_fetch_and_add(&(shared_memory->ready), 1);
#endif

    for (int proc_idx = 0; proc_idx < par.threads; ++proc_idx) {
        int status = 0;
        while (waitpid(-1, &status, 0) == -1) {
            if (errno == EINTR) {
                continue;
            }
        }
    }


    reader.close();

#ifdef HAVE_MPI
    std::pair<std::string, std::string> outDb = Util::createTmpFileNames(par.db2, par.db2Index, MMseqsMPI::rank);
#else
    std::pair<std::string, std::string> outDb = std::make_pair(par.db2, par.db2Index);
#endif

    std::vector<std::pair<std::string, std::string>> splitFiles;
    for (int proc_idx = 0; proc_idx < par.threads; ++proc_idx) {
        splitFiles.emplace_back(Util::createTmpFileNames(outDb.first, outDb.second, proc_idx));
    }
    DBWriter::mergeResults(outDb.first, outDb.second, splitFiles);

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    if (MMseqsMPI::rank == MMseqsMPI::MASTER) {
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
#endif
