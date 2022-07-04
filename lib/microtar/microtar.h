/**
 * Copyright (c) 2017 rxi
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the MIT license. See `microtar.c` for details.
 */

#ifndef MICROTAR_H
#define MICROTAR_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>

#define MTAR_VERSION "0.1.0"

enum {
  MTAR_ESUCCESS     =  0,
  MTAR_EFAILURE     = -1,
  MTAR_EOPENFAIL    = -2,
  MTAR_EREADFAIL    = -3,
  MTAR_ESEEKFAIL    = -5,
  MTAR_EBADCHKSUM   = -6,
  MTAR_ENULLRECORD  = -7,
  MTAR_ENOTFOUND    = -8
};

enum {
  MTAR_TREG   = '0',
  MTAR_TLNK   = '1',
  MTAR_TSYM   = '2',
  MTAR_TCHR   = '3',
  MTAR_TBLK   = '4',
  MTAR_TDIR   = '5',
  MTAR_TFIFO  = '6',
  MTAR_TCONT  = '7',
  MTAR_TGNU_LONGNAME = 'L',
  MTAR_TGNU_LONGLINK = 'K',
  MTAR_TOLDREG = '\0'
};

typedef struct {
  unsigned mode;
  unsigned owner;
  unsigned size;
  unsigned mtime;
  unsigned type;
  char name[100];
  char linkname[100];
} mtar_header_t;


typedef struct mtar_t mtar_t;

struct mtar_t {
  int (*read)(mtar_t *tar, void *data, size_t size);
  int (*seek)(mtar_t *tar, long pos, int whence);
  int (*close)(mtar_t *tar);
  void *stream;
  size_t curr_size;
  int isFinished;
};


const char* mtar_strerror(int err);

int mtar_open(mtar_t *tar, const char *filename);
int mtar_close(mtar_t *tar);

int mtar_read_header(mtar_t *tar, mtar_header_t *h);
int mtar_read_data(mtar_t *tar, void *ptr, size_t size);
int mtar_skip_data(mtar_t *tar);

#ifdef __cplusplus
}
#endif

#endif
