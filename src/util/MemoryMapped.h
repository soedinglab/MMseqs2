// //////////////////////////////////////////////////////////
// MemoryMapped.h
// Copyright (c) 2013 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
//

#pragma once

// define fixed size integer types
#ifdef _MSC_VER
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

#include <string>


/// Portable read-only memory mapping (Windows and Linux)
/** Filesize limited by size_t, usually 2^32 or 2^64 */
class MemoryMapped
{
public:
  /// tweak performance
  enum CacheHint
  {
    Normal,         ///< good overall performance
    SequentialScan, ///< read file only once with few seeks
    RandomAccess    ///< jump around
  };

  /// how much should be mappend
  enum MapRange
  {
    WholeFile = 0   ///< everything ... be careful when file is larger than memory
  };

  /// do nothing, must use open()
  MemoryMapped();
  /// open file, mappedBytes = 0 maps the whole file
  MemoryMapped(const std::string& filename, size_t mappedBytes = WholeFile, CacheHint hint = Normal);
  /// close file (see close() )
  ~MemoryMapped();

  /// open file, mappedBytes = 0 maps the whole file
  bool open(const std::string& filename, size_t mappedBytes = WholeFile, CacheHint hint = Normal);
  /// close file
  void close();

  /// access position, no range checking (faster)
  unsigned char operator[](size_t offset) const;
  /// access position, including range checking
  unsigned char at        (size_t offset) const;

  /// raw access
  const unsigned char* getData() const;

  /// true, if file successfully opened
  bool isValid() const;

  /// get file size
  uint64_t size() const;
  /// get number of actually mapped bytes
  size_t   mappedSize() const;

  /// replace mapping by a new one of the same file, offset MUST be a multiple of the page size
  bool remap(uint64_t offset, size_t mappedBytes);

private:
  /// don't copy object
  MemoryMapped(const MemoryMapped&);
  /// don't copy object
  MemoryMapped& operator=(const MemoryMapped&);

  /// get OS page size (for remap)
  static int getpagesize();

  /// file name
  std::string _filename;
  /// file size
  uint64_t    _filesize;
  /// caching strategy
  CacheHint   _hint;
  /// mapped size
  size_t      _mappedBytes;

  /// define handle
#ifdef _MSC_VER
  typedef void* FileHandle;
  /// Windows handle to memory mapping of _file
  void*       _mappedFile;
#else
  typedef int   FileHandle;
#endif

  /// file handle
  FileHandle  _file;
  /// pointer to the file contents mapped into memory
  void*       _mappedView;
};
