#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <unistd.h>
#ifdef __APPLE__
#include <mach/mach.h>
#endif

// availability estimate: free plus reclaimable inactive pages
inline size_t marvAvailableHostMemory(size_t total) {
#ifdef __APPLE__
    vm_statistics64_data_t vm{};
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
    const mach_port_t host = mach_host_self();
    const kern_return_t result = host_statistics64(host, HOST_VM_INFO64,
        reinterpret_cast<host_info64_t>(&vm), &count);
    mach_port_deallocate(mach_task_self(), host);
    if (result == KERN_SUCCESS)
        return std::min(total, (size_t(vm.free_count) + vm.inactive_count) * size_t(getpagesize()));
#elif defined(__linux__)
    if (FILE* file = std::fopen("/proc/meminfo", "r")) {
        char line[256];
        unsigned long long kb = 0;
        while (std::fgets(line, sizeof(line), file)) {
            if (std::sscanf(line, "MemAvailable: %llu kB", &kb) == 1) {
                std::fclose(file);
                return std::min(total, size_t(kb) * 1024);
            }
        }
        std::fclose(file);
    }
#endif
    return 0;
}
