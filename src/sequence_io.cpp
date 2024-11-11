/*************************************************************************//**
 *
 * @file  single and pairwise reading of sequences from FASTA/FASTQ files
 *
 * (c) 2017-2023 André Müller (mail@andremueller-online.de)
 * MIT License
 *
 *****************************************************************************/

#include <sstream>
#include <limits>
#include <iostream>

#include "sequence_io.h"

using std::string;


//-----------------------------------------------------------------------------
// BATCHED READING 
//-----------------------------------------------------------------------------
sequence_batch
read_all_sequences_from_file (const string& filename, int align) 
{
    sequence_batch batch;

    read_all_sequences_from_file(filename, batch, align);

    return batch;
}


//-------------------------------------------------------------------
void read_all_sequences_from_file (const string& filename, sequence_batch& out, 
                                   int alignment) 
{
    auto reader = make_sequence_reader(filename);

    out.chars.clear();
    out.offsets.clear();
    out.lengths.clear();
    out.headers.clear();
    out.qualities.clear();

    if (!reader) return;

    sequence_reader::data_type buffer;
    while (reader->has_next()) {
        reader->next_data(buffer);
        if (!buffer.empty()) {
            // pad if alignment criterion not met
            const auto miss = out.chars.size() % alignment;
            if (miss > 0) out.chars.insert(out.chars.end(), alignment-miss, ' '); 
            out.offsets.push_back(out.chars.size());
            out.lengths.push_back(buffer.size());
            out.chars.insert(out.chars.end(), buffer.begin(), buffer.end());
        }
    }
    out.offsets.push_back(out.chars.size());
}




//-------------------------------------------------------------------
sequence_batch
read_all_sequences_and_headers_from_file (const string& filename, int align) 
{
    sequence_batch batch;

    read_all_sequences_and_headers_from_file(filename, batch, align);

    return batch;
}


//-------------------------------------------------------------------
void read_all_sequences_and_headers_from_file (const string& filename,
                                               sequence_batch& out, 
                                               int alignment)
{
    auto reader = make_sequence_reader(filename);

    out.chars.clear();
    out.offsets.clear();
    out.lengths.clear();
    out.headers.clear();
    out.qualities.clear();

    if (!reader) return;

    sequence_reader::header_type headerBuf;
    sequence_reader::data_type dataBuf;

    while (reader->has_next()) {
        reader->next_header_and_data(headerBuf, dataBuf);
        if (!dataBuf.empty()) {
            // pad if alignment criterion not met
            const auto miss = out.chars.size() % alignment;
            if (miss > 0) out.chars.insert(out.chars.end(), alignment-miss, ' '); 
            out.offsets.push_back(out.chars.size());
            out.lengths.push_back(dataBuf.size());
            out.chars.insert(out.chars.end(), dataBuf.begin(), dataBuf.end());
            out.headers.push_back(headerBuf);
        }
    }
    out.offsets.push_back(out.chars.size());
}




//-------------------------------------------------------------------
sequence_batch
read_all_sequences_and_meta_info_from_file (const string& filename, int align) 
{
    sequence_batch batch;

    read_all_sequences_and_meta_info_from_file(filename, batch, align);

    return batch;
}


//-------------------------------------------------------------------
void read_all_sequences_and_meta_info_from_file (const string& filename,
                                               sequence_batch& out, 
                                               int alignment)
{
    auto reader = make_sequence_reader(filename);

    out.chars.clear();
    out.offsets.clear();
    out.lengths.clear();
    out.headers.clear();
    out.qualities.clear();

    if (!reader) return;

    sequence_reader::header_type headerBuf;
    sequence_reader::data_type dataBuf;
    sequence_reader::qualities_type qualBuf;

    while (reader->has_next()) {
        reader->next_header_data_qualities(headerBuf, dataBuf, qualBuf);
        if (!dataBuf.empty()) {
            // pad if alignment criterion not met
            const auto miss = out.chars.size() % alignment;
            if (miss > 0) out.chars.insert(out.chars.end(), alignment-miss, ' '); 
            out.offsets.push_back(out.chars.size());
            out.lengths.push_back(dataBuf.size());
            out.chars.insert(out.chars.end(), dataBuf.begin(), dataBuf.end());
            out.headers.push_back(headerBuf);
            out.qualities.push_back(qualBuf);
        }
    }
    out.offsets.push_back(out.chars.size());
}




//-----------------------------------------------------------------------------
// SEQUENCE_READER  B A S E 
//-----------------------------------------------------------------------------
sequence_reader::sequence
sequence_reader::next ()
{
    sequence seq;
    next(seq);
    return seq;
}



//-------------------------------------------------------------------
void sequence_reader::next (sequence& seq)
{
    if (!has_next()) return;

    ++index_;
    seq.index = index_;
    read_next(&seq.header, &seq.data, &seq.qualities);
}



//-------------------------------------------------------------------
sequence_reader::header_type
sequence_reader::next_header ()
{
    if (!has_next()) return header_type{};

    ++index_;
    header_type header;
    read_next(&header, nullptr, nullptr);
    return header;
}



//-------------------------------------------------------------------
sequence_reader::data_type
sequence_reader::next_data ()
{
    if (!has_next()) return data_type{};

    ++index_;
    data_type data;
    read_next(nullptr, &data, nullptr);
    return data;
}



//-------------------------------------------------------------------
sequence_reader::index_type
sequence_reader::next_data (sequence::data_type& data)
{
    if (!has_next()) {
        data.clear();
        return index(); 
    }

    ++index_;
    read_next(nullptr, &data, nullptr);
    return index_;
}



//-------------------------------------------------------------------
sequence_reader::index_type
sequence_reader::next_header_and_data (sequence::header_type& header,
                                       sequence::data_type& data)
{
    if (!has_next()) {
        header.clear();
        data.clear();
        return index();
    }

    ++index_;
    read_next(&header, &data, nullptr);
    return index_;
}



//-------------------------------------------------------------------
sequence_reader::index_type
sequence_reader::next_data_and_qualities (
    sequence::data_type& data, sequence::qualities_type& qual)
{
    if (!has_next()) {
        data.clear();
        qual.clear();
        return index();
    }

    ++index_;
    read_next(nullptr, &data, &qual);
    return index_;
}



//-------------------------------------------------------------------
sequence_reader::index_type
sequence_reader::next_header_data_qualities (
    sequence::header_type& header,
    sequence::data_type& data,
    sequence::qualities_type& qual)
{
    if (!has_next()) {
        header.clear();
        data.clear();
        qual.clear();
        return index();
    }

    ++index_;
    read_next(&header, &data, &qual);
    return index_;
}



//-------------------------------------------------------------------
void sequence_reader::skip (index_type skip)
{
    if (skip < 1) return;

    for(; skip > 0 && has_next(); --skip) {
        ++index_;
        skip_next();
    }
}






//-----------------------------------------------------------------------------
// F A S T A    R E A D E R
//-----------------------------------------------------------------------------
fasta_reader::fasta_reader (const string& filename):
    sequence_reader{},
    file_{},
    linebuffer_{},
    pos_{0}
{
    if (!filename.empty()) {
        file_.open(filename);

        if (!file_.good()) {
            invalidate();
            throw file_access_error{"can't open file " + filename};
        }
    }
    else {
        throw file_access_error{"no filename was given"};
    }
}



//-------------------------------------------------------------------
void fasta_reader::read_next (header_type* header, data_type* data, qualities_type*)
{
    if (linebuffer_.empty()) {
        getline(file_, linebuffer_);
    }
    pos_ += linebuffer_.size() + 1;

    if (linebuffer_[0] != '>') {
        throw io_format_error{"malformed fasta file - expected header char > not found"};
        invalidate();
        return;
    }

    if (header) *header = linebuffer_.substr(1);

    if (data) data->clear();

    while (file_.good()) {
        getline(file_, linebuffer_);
        if (linebuffer_[0] == '>') {
            break;
        }
        else {
            if (data) data->append(linebuffer_);
            pos_ += linebuffer_.size() + 1;
        }
    }

    if (data && data->empty()) {
        throw io_format_error{"malformed fasta file - zero-length sequence"
                              + (header ? *header : header_type{""})};
        invalidate();
        return;
    }

    if (!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}




//-------------------------------------------------------------------
void fasta_reader::skip_next ()
{
    if (linebuffer_.empty()) {
        file_.ignore(1);
        pos_ += file_.gcount();
    } else {
        pos_ += linebuffer_.size() + 1;
        linebuffer_.clear();
    }
    file_.ignore(std::numeric_limits<std::streamsize>::max(), '>');
    pos_ += file_.gcount();

    if (file_.good()) {
        file_.unget();
        pos_ -= 1;
    } else {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
void fasta_reader::do_seek (std::streampos pos)
{
    file_.seekg(pos);
    pos_ = pos;
    linebuffer_.clear();

    if (!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
std::streampos fasta_reader::do_tell ()
{
    return pos_;
}






//-----------------------------------------------------------------------------
// F A S T Q    R E A D E R
//-----------------------------------------------------------------------------
fastq_reader::fastq_reader (const string& filename):
    sequence_reader{},
    file_{}, linebuffer_{}, pos_{0}
{
    if (!filename.empty()) {
        file_.open(filename);

        if (!file_.good()) {
            invalidate();
            throw file_access_error{"can't open file " + filename};
        }
    }
    else {
        throw file_access_error{"no filename was given"};
    }
}



//-------------------------------------------------------------------
void fastq_reader::read_next (header_type* header, data_type* data,
                              qualities_type* qualities)
{
    // 1st line (data header)
    getline(file_, linebuffer_);

    if (linebuffer_.empty()) {
        pos_ = -1;
        invalidate();
        if (header) header->clear();
        if (data) data->clear();
        if (qualities) qualities->clear();
        return;
    }

    pos_ += linebuffer_.size() + 1;

    if (linebuffer_[0] != '@') {
        if (linebuffer_[0] != '\r') {
            throw io_format_error{"malformed fastq file - sequence header: "  + linebuffer_};
        }
        invalidate();
        return;
    }

    if (header) {
        *header = linebuffer_.substr(1);
    }

    // 2nd line (sequence data)
    if (data) {
        getline(file_, *data);
        pos_ += data->size() + 1;
    } else {
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        pos_ += file_.gcount();
    }

    // 3rd (qualities header) + 4th line (qualities)
    if (qualities) {
        getline(file_, linebuffer_);
        pos_ += linebuffer_.size() + 1;

        if (linebuffer_.empty() || linebuffer_[0] != '+') {
            if (linebuffer_[0] != '\r') {
                throw io_format_error{"malformed fastq file - quality header: "  + linebuffer_};
            }
            invalidate();
            return;
        }

        getline(file_, *qualities);
        pos_ += qualities->size() + 1;
    }
    else {
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        pos_ += file_.gcount();
        file_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        pos_ += file_.gcount();
    }

    if (!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
void fastq_reader::skip_next ()
{
    read_next(nullptr, nullptr, nullptr);
}



//-------------------------------------------------------------------
void fastq_reader::do_seek (std::streampos pos)
{
    file_.seekg(pos);
    pos_ = pos;

    if (!file_.good()) {
        pos_ = -1;
        invalidate();
    }
}



//-------------------------------------------------------------------
std::streampos fastq_reader::do_tell ()
{
    return pos_;
}






//-----------------------------------------------------------------------------
// P A I R    R E A D E R
//-----------------------------------------------------------------------------
sequence_pair_reader::sequence_pair_reader (const std::string& filename1,
                                            const std::string& filename2)
:
    reader1_{nullptr},
    reader2_{nullptr},
    singleMode_{true}
{
    if (!filename1.empty()) {
        reader1_ = make_sequence_reader(filename1);

        if (!filename2.empty()) {
            singleMode_ = false;
            if (filename1 != filename2) {
                reader2_ = make_sequence_reader(filename2);
            }
        }
    }
}



//-------------------------------------------------------------------
bool sequence_pair_reader::has_next () const noexcept
{
    if (!reader1_) return false;
    if (!reader1_->has_next()) return false;
    if (!reader2_) return true;
    if (!reader2_->has_next()) return false;
    return true;
}



//-------------------------------------------------------------------
sequence_pair_reader::sequence_pair
sequence_pair_reader::next ()
{
    sequence_pair seq;
    next(seq);
    return seq;
}



//-------------------------------------------------------------------
void sequence_pair_reader::next (sequence_pair& seq)
{
    if (!has_next()) return;

    // only one sequence per call
    if (singleMode_) {
        reader1_->next(seq.first);
        seq.second.header.clear();
        seq.second.data.clear();
        seq.second.qualities.clear();
    }
    // pair = single sequences from 2 separate files (read in lockstep)
    else if (reader2_) {
        reader1_->next(seq.first);
        reader2_->next(seq.second);
    }
    // pair = 2 consecutive sequences from same file
    else {
        const auto idx = reader1_->index();
        reader1_->next(seq.first);
        //make sure the index is only increased after the 2nd 'next()'
        reader1_->index_offset(idx);
        reader1_->next(seq.second);
    }
}



//-------------------------------------------------------------------
sequence_pair_reader::header_type
sequence_pair_reader::next_header ()
{
    if (!has_next()) return header_type{};

    // only one sequence per call
    if (singleMode_) {
        return reader1_->next_header();
    }

    // pair = single sequences from 2 separate files (read in lockstep)
    if (reader2_) {
        reader2_->next_header();
        return reader1_->next_header();
    }

    // pair = 2 consecutive sequences from same file
    const auto idx = reader1_->index();
    auto header = reader1_->next_header();
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    reader1_->next_header();
    return header;
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type
sequence_pair_reader::next_data (sequence::data_type& data1,
                                 sequence::data_type& data2)
{
    if (!has_next()) return index();

    // only one sequence per call
    if (singleMode_) {
        data2.clear();
        return reader1_->next_data(data1);
    }

    // pair = single sequences from 2 separate files (read in lockstep)
    if (reader2_) {
        reader1_->next_data(data1);
        return reader2_->next_data(data2);
    }

    // pair = 2 consecutive sequences from same file
    const auto idx = reader1_->index();
    reader1_->next_data(data1);
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    return reader1_->next_data(data2);
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type
sequence_pair_reader::next_header_and_data (sequence::header_type& header1,
                                            sequence::data_type& data1,
                                            sequence::data_type& data2)
{
    if (!has_next()) return index();

    // only one sequence per call
    if (singleMode_) {
        data2.clear();
        return reader1_->next_header_and_data(header1, data1);
    }

    // pair = single sequences from 2 separate files (read in lockstep)
    if (reader2_) {
        reader1_->next_header_and_data(header1, data1);
        return reader2_->next_data(data2);
    }

    // pair = 2 consecutive sequences from same file
    const auto idx = reader1_->index();
    reader1_->next_header_and_data(header1, data1);
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    return reader1_->next_data(data2);
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type
sequence_pair_reader::next_data_and_qualities (
    sequence::data_type& data1,
    sequence::data_type& data2,
    sequence::qualities_type& qual1,
    sequence::qualities_type& qual2)
{
    if (!has_next()) return index();

    // only one sequence per call
    if (singleMode_) {
        data2.clear();
        return reader1_->next_data_and_qualities(data1, qual1);
    }

    // pair = single sequences from 2 separate files (read in lockstep)
    if (reader2_) {
        reader1_->next_data_and_qualities(data1, qual1);
        return reader2_->next_data_and_qualities(data2, qual2);
    }

    // pair = 2 consecutive sequences from same file
    const auto idx = reader1_->index();
    reader1_->next_data_and_qualities(data1, qual1);
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    return reader1_->next_data_and_qualities(data2, qual2);
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type
sequence_pair_reader::next_header_data_qualities (
    sequence::header_type& header1,
    sequence::data_type& data1,
    sequence::data_type& data2,
    sequence::qualities_type& qual1,
    sequence::qualities_type& qual2)
{
    if (!has_next()) return index();

    // only one sequence per call
    if (singleMode_) {
        data2.clear();
        return reader1_->next_header_data_qualities(header1, data1, qual1);
    }

    // pair = single sequences from 2 separate files (read in lockstep)
    if (reader2_) {
        reader1_->next_header_data_qualities(header1, data1, qual1);
        return reader2_->next_data_and_qualities(data2, qual2);
    }

    // pair = 2 consecutive sequences from same file
    const auto idx = reader1_->index();
    reader1_->next_header_data_qualities(header1, data1, qual1);
    //make sure the index is only increased after the 2nd 'next()'
    reader1_->index_offset(idx);
    return reader1_->next_data_and_qualities(data2, qual2);
}



//-------------------------------------------------------------------
void sequence_pair_reader::skip (index_type skip)
{
    if (skip < 1 || !reader1_) return;

    if (reader2_) {
        reader1_->skip(skip);
        reader2_->skip(skip);
    }
    else if (singleMode_) {
        reader1_->skip(skip);
    }
    else {
        const auto idx = reader1_->index();
        reader1_->skip(2*skip);
        reader1_->index_offset(idx+skip);
    }
}



//-------------------------------------------------------------------
sequence_pair_reader::index_type sequence_pair_reader::index () const noexcept
{
    if (!reader1_) return index_type{0};
    return reader1_->index();
}



//-------------------------------------------------------------------
void sequence_pair_reader::index_offset (index_type index)
{
    if (!reader1_) return;

    reader1_->index_offset(index);
    if (reader2_) reader2_->index_offset(index);
}



//-------------------------------------------------------------------
void sequence_pair_reader::seek (const stream_positions& pos)
{
    if (!reader1_) return;

    reader1_->seek(pos.first);
    if (reader2_) reader2_->seek(pos.second);
}



//-------------------------------------------------------------------
sequence_pair_reader::stream_positions
sequence_pair_reader::tell ()
{
    return stream_positions{
            reader1_ ? reader1_->tell() : std::streampos{},
            reader2_ ? reader2_->tell() : std::streampos{} };
}






//-------------------------------------------------------------------
std::unique_ptr<sequence_reader>
make_sequence_reader (const string& filename)
{
    if (filename.empty()) return nullptr;

    auto n = filename.size();
    if (filename.find(".fq")    == (n-3) ||
       filename.find(".fnq")   == (n-4) ||
       filename.find(".fastq") == (n-6) )
    {
        return std::make_unique<fastq_reader>(filename);
    }
    else if (filename.find(".fa")    == (n-3) ||
            filename.find(".fna")   == (n-4) ||
            filename.find(".fasta") == (n-6) )
    {
        return std::make_unique<fasta_reader>(filename);
    }

    //try to determine file type content
    std::ifstream is {filename};
    if (is.good()) {
        string line;
        getline(is,line);
        if (!line.empty()) {
            if (line[0] == '>') {
                return std::make_unique<fasta_reader>(filename);
            }
            else if (line[0] == '@') {
                return std::make_unique<fastq_reader>(filename);
            }
        }
        throw file_read_error{"file format not recognized"};
    }

    throw file_access_error{"file not accessible"};
    return nullptr;
}

