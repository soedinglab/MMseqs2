/* Adapted from  MemoryMapped.h (http://create.stephan-brumme.com/portable-memory-mapping/)
 * Copyright (c) 2013 Stephan Brumme. All rights reserved.
 * License
 * All source code published on http://create.stephan-brumme.com and its sub-pages is licensed similar to the zlib license:
 *
 * This software is provided 'as-is', without any express or implied warranty. In no event will the author be held liable for any damages arising from the use of this software.
 * Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:
 * The origin of this software must not be misrepresented; you must not claim that you wrote the original software.
 * If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
 * Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
*/

#include "MemoryMapped.h"

#include <cstdio>
#include "Debug.h"
#include "Util.h"

// OS-specific
#ifdef _MSC_VER
// Windows
#include <windows.h>

#define MM_INVALID_FILE NULL
#else
// Linux
// enable large file support on 32 bit systems
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#ifdef  _FILE_OFFSET_BITS
#undef  _FILE_OFFSET_BITS
#endif
#define _FILE_OFFSET_BITS 64
// and include needed headers
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>

#define MM_INVALID_FILE -1
#endif


/// do nothing, must use open()
MemoryMapped::MemoryMapped()
        : _filename   (),
          _filesize   (0),
          _hint       (Normal),
          _mappedBytes(0),
          _file       (MM_INVALID_FILE),
#ifdef _MSC_VER
        _mappedFile (NULL),
#endif
          _mappedView (NULL),
          _emptyView ("") {}


/// open file, mappedBytes = 0 maps the whole file
MemoryMapped::MemoryMapped(const std::string& filename, size_t mappedBytes, CacheHint hint)
        : _filename   (filename),
          _filesize   (0),
          _hint       (hint),
          _mappedBytes(mappedBytes),
          _file       (MM_INVALID_FILE),
#ifdef _MSC_VER
        _mappedFile (NULL),
#endif
          _mappedView (NULL),
          _emptyView ("") {
    open(filename, mappedBytes, hint);
}


/// close file (see close() )
MemoryMapped::~MemoryMapped()
{
    close();
}


/// open file
bool MemoryMapped::open(const std::string& filename, size_t mappedBytes, CacheHint hint)
{
    // already open ?
    if (isValid())
        return false;

    _file       = MM_INVALID_FILE;
    _filesize   = 0;
    _hint       = hint;
#ifdef _MSC_VER
    _mappedFile = NULL;
#endif
    _mappedView = NULL;

#ifdef _MSC_VER
    // Windows

  DWORD winHint = 0;
  switch (_hint)
  {
  case Normal:         winHint = FILE_ATTRIBUTE_NORMAL;     break;
  case SequentialScan: winHint = FILE_FLAG_SEQUENTIAL_SCAN; break;
  case RandomAccess:   winHint = FILE_FLAG_RANDOM_ACCESS;   break;
  default: break;
  }

  // open file
  _file = ::CreateFileA(filename.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, winHint, NULL);
  if (_file == MM_INVALID_FILE)
    return false;

  // file size
  LARGE_INTEGER result;
  if (!GetFileSizeEx(_file, &result))
    return false;
  _filesize = static_cast<uint64_t>(result.QuadPart);

  // convert to mapped mode
  _mappedFile = ::CreateFileMapping(_file, NULL, PAGE_READONLY, 0, 0, NULL);
  if (!_mappedFile)
    return false;

#else

    // Linux

    // open file

    _file = ::open(filename.c_str(), O_RDONLY);
    if (_file == MM_INVALID_FILE)
    {
        return false;
    }

    // file size
    struct stat statInfo;
    if (fstat(_file, &statInfo) < 0)
        return false;

    _filesize = statInfo.st_size;
#endif

    // initial mapping
    remap(0, mappedBytes);

    if (!_mappedView)
        return false;

    // everything's fine
    return true;
}


/// close file
void MemoryMapped::close()
{
    // kill pointer
    if (_mappedView) {
        if (_mappedView != _emptyView) {
#ifdef _MSC_VER
            ::UnmapViewOfFile(_mappedView);
#else
            ::munmap(_mappedView, _filesize);
#endif
        }
        _mappedView = NULL;
    }

#ifdef _MSC_VER
    if (_mappedFile)
  {
    ::CloseHandle(_mappedFile);
    _mappedFile = NULL;
  }
#endif

    // close underlying file
    if (_file != MM_INVALID_FILE)
    {
#ifdef _MSC_VER
        ::CloseHandle(_file);
#else
        ::close(_file);
#endif
        _file = MM_INVALID_FILE;
    }

    _filesize = 0;
}


/// access position, no range checking (faster)
unsigned char MemoryMapped::operator[](size_t offset) const
{
    return ((unsigned char*)_mappedView)[offset];
}


/// access position, including range checking
unsigned char MemoryMapped::at(size_t offset) const
{
    // checks
    if (offset >= _filesize) {
            Debug(Debug::ERROR) << "offset (" << offset << ") >= _filesize (" << _filesize << ")\n";
            EXIT(EXIT_FAILURE);
    }
    return operator[](offset);
}


/// raw access
const unsigned char* MemoryMapped::getData() const
{
    return (const unsigned char*)_mappedView;
}


/// true, if file successfully opened
bool MemoryMapped::isValid() const
{
    return _mappedView != NULL;
}


/// get file size
uint64_t MemoryMapped::size() const
{
    return _filesize;
}


/// get number of actually mapped bytes
size_t MemoryMapped::mappedSize() const
{
    return _mappedBytes;
}


/// replace mapping by a new one of the same file, offset MUST be a multiple of the page size
bool MemoryMapped::remap(uint64_t offset, size_t mappedBytes)
{
    if (_file == MM_INVALID_FILE)
        return false;

    if (mappedBytes == WholeFile)
        mappedBytes = _filesize;

    // close old mapping
    if (_mappedView) {
        if (_mappedView != _emptyView) {
#ifdef _MSC_VER
            ::UnmapViewOfFile(_mappedView);
#else
            ::munmap(_mappedView, _mappedBytes);
#endif
        }
        _mappedView = NULL;
    }

    // don't go further than end of file
    if (offset > _filesize)
        return false;
    if (offset + mappedBytes > _filesize)
        mappedBytes = size_t(_filesize - offset);

    // mmap doesn't opening 0-byte large files
    // map a fake null view to support this case
    if (mappedBytes == 0) {
        _mappedView = (void*)_emptyView;
        return true;
    }

#ifdef _MSC_VER
    // Windows

  DWORD offsetLow  = DWORD(offset & 0xFFFFFFFF);
  DWORD offsetHigh = DWORD(offset >> 32);
  _mappedBytes = mappedBytes;

  // get memory address
  _mappedView = ::MapViewOfFile(_mappedFile, FILE_MAP_READ, offsetHigh, offsetLow, mappedBytes);

  if (_mappedView == NULL)
  {
    _mappedBytes = 0;
    _mappedView  = NULL;
    return false;
  }

  return true;

#else

    // Linux
    // new mapping
    _mappedView = ::mmap(NULL, mappedBytes, PROT_READ, MAP_SHARED, _file, offset);
    if (_mappedView == MAP_FAILED)
    {
        _mappedBytes = 0;
        _mappedView  = NULL;
        return false;
    }

    _mappedBytes = mappedBytes;

    // tweak performance
    int linuxHint = 0;
    switch (_hint)
    {
        case Normal:         linuxHint = MADV_NORMAL;     break;
        case SequentialScan: linuxHint = MADV_SEQUENTIAL; break;
        case RandomAccess:   linuxHint = MADV_RANDOM;     break;
        default: break;
    }
    // assume that file will be accessed soon
    //linuxHint |= MADV_WILLNEED;
    // assume that file will be large
    //linuxHint |= MADV_HUGEPAGE;

    ::madvise(_mappedView, _mappedBytes, linuxHint);

    return true;
#endif
}


/// get OS page size (for remap)
int MemoryMapped::getpagesize()
{
#ifdef _MSC_VER
    SYSTEM_INFO sysInfo;
  GetSystemInfo(&sysInfo);
  return sysInfo.dwAllocationGranularity;
#else
    return sysconf(_SC_PAGESIZE); //::getpagesize();
#endif
}
