/*
 * Based on https://github.com/TuftsBCB/io
 * Ported from GO to C++
 * License:
 * Copyright (c) <2013> <Andrew Gallant>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */
#include "A3MReader.h"

#include <sstream>
#include <algorithm>
#include <iterator>

#include "kseq.h"
#include "KSeqBufferReader.h"

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

A3mReader::A3mReader(std::string a3m) : length(0) {
    kseq_buffer_t buffer(const_cast<char*>(a3m.c_str()), a3m.length());
    kseq_t *seq = kseq_init(&buffer);
    while (kseq_read(seq) >= 0) {
        std::string header = seq->name.s;
        header.append(" ");
        header.append(seq->comment.s);
        headers.push_back(header);

        std::string sequence = seq->seq.s;
        addSequence(sequence);
    }
    kseq_destroy(seq);
}

std::string A3mReader::getFasta() {
    if(entries.size() < 1) {
        return std::string();
    }

    const std::vector<char>& query = entries[0];

    std::ostringstream ss;
    for (size_t i = 0; i < entries.size(); ++i) {
        ss << ">" << headers[i] << "\n";

        const std::vector<char> &entry = entries[i];
        for (size_t i = 0; i < query.size(); ++i) {
            if(query[i] == '.' || query[i] == '-') {
                continue;
            } else {
                if(entry[i] == '.') {
                    ss << '-';
                } else {
                    ss << entry[i];
                }
            }
        }

        ss << "\n";
    }

    return ss.str();
}

enum HMMState {
    Deletion = 0,
    Insertion,
    Match
};

HMMState getState(char r) {
    if (r >= 'a' && r <= 'z') {
        return Insertion;
    } else if (r >= 'A' && r <= 'Z') {
        return Match;
    }

    switch (r) {
        case '-':
            return Deletion;
        case '.':
            return Insertion;
        default:
            return Match;
    }
}

void addInsert(std::vector<char>& sequence, size_t col) {
    sequence.insert(sequence.begin() + col, '.');
}

char translateA2M(char b) {
    if ((b >= 'a' && b <= 'z') || (b >= 'A' && b <= 'Z')) {
        return b;
    }

    switch (b) {
        case '*':
            return 0;
        case '-':
            return '-';
        case '.':
            return '.';
        case '/':
            return '-';
        default:
            return 0;
    }
}

bool A3mReader::columnHasInsertion(size_t col) {
    for (std::vector<std::vector<char>>::iterator it = entries.begin(); it != entries.end(); ++it) {
        if (getState((*it)[col]) == Insertion) {
            return true;
        }
    }
    return false;
}

void A3mReader::addSequence(const std::string& sequence) {
    if (sequence.empty()) {
        return;
    }


    // We will, in all likelihood, modify the sequence if it's in A3M format.
    // So we copy it to prevent weird effects to the caller.
    std::vector<char> copy(sequence.begin(), sequence.end());
    for (std::vector<char>::iterator it = copy.begin(); it != copy.end(); ++it) {
        (*it) = translateA2M(*it);
    }

    // The first sequence is easy.
    if (entries.empty()) {
        entries.push_back(copy);
        length = sequence.size();
        return;
    }

    // This should be an A3M formatted sequence (no way to do a sanity check
    // though, I don't think).
    // Therefore, we need to assimilate it with other sequences, which will
    // change the length of the MSA. In particular, we'll need to add '.'
    // residues to existing sequences.
    for (size_t col = 0; col < length; ++col) {
        bool colHasInsert = columnHasInsertion(col);
        if (col >= copy.size()) {
            if (colHasInsert) {
                copy.push_back('.');
            } else {
                copy.push_back('-');
            }
            continue;
        }

        bool seqHasInsert = getState(copy[col]) == Insertion;
        if (colHasInsert && seqHasInsert) {
            // do nothing, we're in sync
        } else if (colHasInsert) {
            // Put an insert into the sequence we're adding.
            addInsert(copy, col);
        } else if (seqHasInsert) {
            // Put an insert into the rest of the sequences.
            for (std::vector<std::vector<char>>::iterator it = entries.begin(); it != entries.end(); ++it) {
                addInsert(*it, col);
            }

            length++;
        }
    }
    entries.push_back(copy);
}
