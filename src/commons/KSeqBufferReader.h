#ifndef KSEQ_BUFFER_READER_H
#define KSEQ_BUFFER_READER_H

#include <sys/types.h>

typedef struct kseq_buffer {
    char* buffer;
    size_t position, length;

    kseq_buffer () : buffer(NULL), position(0), length(0) {};
    kseq_buffer (char* buffer, size_t length) : buffer(buffer), position(0), length(length) {};
} kseq_buffer_t;

inline ssize_t kseq_buffer_reader(kseq_buffer_t *inBuffer, char *outBuffer, size_t nbyte) {
    if (inBuffer->position > inBuffer->length) {
        return 0;
    }

    size_t bytes = nbyte;
    if (inBuffer->position + bytes > inBuffer->length) {
        bytes = inBuffer->length - inBuffer->position;
    }

    if (bytes == 0) {
        return 0;
    }

    for (size_t i = inBuffer->position; i < inBuffer->position + bytes; ++i) {
        size_t index = i - inBuffer->position;
        outBuffer[index] = inBuffer->buffer[i];
    }

    inBuffer->position += bytes;

    return bytes;
}

#endif
