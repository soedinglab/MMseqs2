/* The MIT License

   Copyright (c) Milot Mirdita

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#ifndef KSEQ_BUFFER_READER_H
#define KSEQ_BUFFER_READER_H

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
