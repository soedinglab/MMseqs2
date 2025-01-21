/*************************************************************************//**
 *
 * @file  single and pairwise reading of sequences from FASTA/FASTQ files
 *
 * (c) 2017-2023 André Müller (mail@andremueller-online.de)
 * MIT License
 *
 *****************************************************************************/

#ifndef SEQUENCE_IO_H_
#define SEQUENCE_IO_H_


#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>



/*************************************************************************//**
 *
 * @member  chars    concatenation of all sequences
 * @member  offsets  start indices of individual sequences
 *
 *****************************************************************************/
struct sequence_batch {
    std::vector<char> chars;               
    std::vector<std::size_t> offsets;  
    std::vector<std::size_t> lengths;  
    std::vector<std::string> headers;  
    std::vector<std::string> qualities;  
};






/*************************************************************************//**
 *
 * @brief  reads ALL sequences from FASTA/FASTQ file, discards headers/qualities;
 *         allocates new sequence_batch
 *
 * @param  align:  alignment of sequence characters in bytes
 *
 *****************************************************************************/
sequence_batch
read_all_sequences_from_file (const std::string& filename, 
                              int charAlignment = 4);



/*************************************************************************//**
 *
 * @brief  reads ALL sequences from FASTA/FASTQ file, discards headers/qualities;
 *         (re-)uses external sequence_batch buffer
 *
 * @param  align:  alignment of sequence characters in bytes
 *
 *****************************************************************************/
void read_all_sequences_from_file (const std::string& filename,
                                   sequence_batch&,
                                   int charAlignment = 4);






/*************************************************************************//**
 *
 * @brief  reads ALL sequences and headers from FASTA/FASTQ file,
 *         discards FASTQ qualities;
 *         allocates new sequence_batch
 *
 * @param  align:  alignment of sequence characters in bytes
 *
 *****************************************************************************/
sequence_batch 
read_all_sequences_and_headers_from_file (const std::string& filename,
                                          int charAlignment = 4);



/*************************************************************************//**
 *
 * @brief  reads ALL sequences from FASTA/FASTQ file, discards FASTQ qualities;
 *         (re-)uses external sequence_batch buffer
 *
 * @param  align:  alignment of sequence characters in bytes
 *
 *****************************************************************************/
void read_all_sequences_and_headers_from_file (const std::string& filename, 
                                               sequence_batch&,
                                               int charAlignment = 4);






/*************************************************************************//**
 *
 * @brief  reads ALL sequences and headers from FASTA/FASTQ file,
 *         allocates new sequence_batch
 *
 * @param  align:  alignment of sequence characters in bytes
 *
 *****************************************************************************/
sequence_batch 
read_all_sequences_and_meta_info_from_file (const std::string& filename,
                                            int charAlignment = 4);



/*************************************************************************//**
 *
 * @brief  reads ALL sequences from FASTA/FASTQ file
 *         (re-)uses external sequence_batch buffer
 *
 * @param  align:  alignment of sequence characters in bytes
 *
 *****************************************************************************/
void read_all_sequences_and_meta_info_from_file (const std::string& filename, 
                                                 sequence_batch&,
                                                 int charAlignment = 4);




/*************************************************************************//**
 *
 * exception types
 *
 *****************************************************************************/
struct io_error : public std::runtime_error {
    explicit io_error (const std::string& what): std::runtime_error(what) {}
};


struct io_format_error : public io_error {
    explicit io_format_error (const std::string& what): io_error(what) {}
};


struct file_io_error : public io_error {
    explicit
    file_io_error (const std::string& what): io_error(what) {}

    explicit
    file_io_error(const std::string& what, const std::string& filename):
        io_error(what), filename_(filename)
    {}

    const char* filename() const noexcept { return filename_.c_str(); }

private:
    std::string filename_;
};


struct file_access_error : public file_io_error {
    explicit 
    file_access_error (const std::string& what): file_io_error(what) {}

    explicit
    file_access_error (const std::string& what, const std::string& filename):
        file_io_error(what,filename)
    {}
};

struct file_read_error : public file_io_error {
    explicit file_read_error (const std::string& what): file_io_error(what) {}

    explicit
    file_read_error (const std::string& what, const std::string& filename):
        file_io_error(what,filename)
    {}
};




/*************************************************************************//**
 *
 * @brief polymorphic file reader for (pairs of) bio-sequences
 *        NOT concurrency safe
 *
 *****************************************************************************/
class sequence_reader
{
public:
    using index_type     = std::uint_least64_t;
    using header_type    = std::string;
    using data_type      = std::string;
    using qualities_type = std::string;

    /** @brief data associated with one sequence */
    struct sequence {
        using index_type     = sequence_reader::index_type;
        using header_type    = sequence_reader::header_type;
        using data_type      = sequence_reader::data_type;
        using qualities_type = sequence_reader::qualities_type;

        index_type  index;      // number of sequence in file (+ offset)
        header_type header;     // meta information (FASTA >, FASTQ @)
        data_type   data;       // actual sequence data
        data_type   qualities;  // quality scores (FASTQ)
    };


    sequence_reader(): index_{0}, valid_{true} {}

    sequence_reader (const sequence_reader&) = delete;
    sequence_reader& operator = (const sequence_reader&) = delete;
    sequence_reader& operator = (sequence_reader&&) = delete;

    virtual ~sequence_reader () = default;


    /** @brief read & return next sequence */
    sequence next ();

    /** @brief read next header only */
    header_type next_header( );

    /** @brief read next sequence data only */
    data_type next_data ();


    /** @brief read next sequence re-using external storage */
    void next (sequence&);

    /** @brief read next header only, re-uses external storage */
    index_type next_header (header_type&);

    /** @brief read next sequence data only, re-uses external storage */
    index_type next_data (data_type&);

    /** @brief read next sequence data & header, re-uses external storage */
    index_type next_header_and_data (header_type&, data_type&);
    
    /** @brief read next header from 1st sequence and data from both sequences
               re-using external storage */
    index_type next_data_and_qualities (sequence::data_type&,
                                        sequence::qualities_type&);

    /** @brief read next sequence data, header & qualities, re-uses external storage */
    index_type next_header_data_qualities (header_type&, data_type&, qualities_type&);


    /** @brief skip n sequences */
    void skip (index_type n);

    bool has_next () const noexcept { return valid_; }

    index_type index () const noexcept { return index_; }

    void index_offset (index_type index) { index_ = index; }

    void seek (std::streampos pos) { do_seek(pos); }
    std::streampos tell ()         { return do_tell(); }


protected:
    void invalidate () { valid_ = false; }


private:
    // derived readers have to implement these
    virtual std::streampos do_tell () = 0;

    virtual void do_seek (std::streampos) = 0;

    virtual void read_next (header_type*, data_type*, qualities_type*) = 0;

    virtual void skip_next () = 0;

    index_type index_;
    bool valid_;
};




/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
class fasta_reader :
    public sequence_reader
{
public:
    explicit
    fasta_reader (const std::string& filename);

private:
    std::streampos do_tell () override;

    void do_seek (std::streampos) override;
    void read_next (header_type*, data_type*, qualities_type*) override;
    void skip_next () override;

private:
    std::ifstream file_;
    std::string linebuffer_;
    std::streampos pos_;
};




/*************************************************************************//**
 *
 *
 *
 *****************************************************************************/
class fastq_reader :
    public sequence_reader
{
public:
    explicit
    fastq_reader (const std::string& filename);

private:
    std::streampos do_tell () override;

    void do_seek (std::streampos) override;
    void read_next (header_type*, data_type*, qualities_type*) override;
    void skip_next () override;

private:
    std::ifstream file_;
    std::string linebuffer_;
    std::streampos pos_;
};





/*************************************************************************//**
 *
 * @brief file reader for (pairs of) bio-sequences
 *        NOT concurrency safe
 *
 *****************************************************************************/
class sequence_pair_reader
{
public:
    using index_type       = sequence_reader::index_type;
    using data_type        = sequence_reader::data_type;
    using header_type      = sequence_reader::header_type;
    using sequence         = sequence_reader::sequence;

    using stream_positions = std::pair<std::streampos,std::streampos>;
    using sequence_pair    = std::pair<sequence,sequence>;


    /** @brief if filename2 empty : single sequence mode
     *         if filename1 == filename2 : read consecutive pairs in one file
     *         else : read from 2 files in lockstep
     */
    sequence_pair_reader (const std::string& filename1,
                          const std::string& filename2);

    sequence_pair_reader (const sequence_pair_reader&) = delete;
    sequence_pair_reader& operator = (const sequence_pair_reader&) = delete;
    sequence_pair_reader& operator = (sequence_pair_reader&&) = delete;

    ~sequence_pair_reader () = default;


    /** @brief read & return next sequence */
    sequence_pair next ();

    /** @brief read next header only */
    header_type next_header ();

    /** @brief read next sequence re-using external storage */
    void next (sequence_pair&);

    /** @brief read next header only, re-uses external storage */
    index_type next_header (sequence::header_type&);

    /** @brief read next sequence data only, re-uses external storage */
    index_type next_data (sequence::data_type&, sequence::data_type&);

    /** @brief read next header from 1st sequence and data from both sequences
               re-using external storage */
    index_type next_header_and_data (sequence::header_type&,
                                     sequence::data_type&,
                                     sequence::data_type&);

    /** @brief read next header from 1st sequence and data from both sequences
               re-using external storage */
    index_type next_data_and_qualities (sequence::data_type&,
                                        sequence::data_type&,
                                        sequence::qualities_type&,
                                        sequence::qualities_type&);

    /** @brief read next header from 1st sequence and 
               data and qualities from both sequences re-using external storage */
    index_type next_header_data_qualities (sequence::header_type&,
                                           sequence::data_type&,
                                           sequence::data_type&,
                                           sequence::qualities_type&,
                                           sequence::qualities_type&);


    /** @brief skip n sequences */
    void skip (index_type n);

    bool has_next () const noexcept;

    index_type index () const noexcept;

    void index_offset (index_type index);

    void seek (const stream_positions& pos);
    stream_positions tell ();


private:
    std::unique_ptr<sequence_reader> reader1_;
    std::unique_ptr<sequence_reader> reader2_;
    bool singleMode_;
};



/*************************************************************************//**
 *
 * @brief guesses and returns a suitable sequence reader
 *        based on a filename pattern
 *
 *****************************************************************************/
std::unique_ptr<sequence_reader>
make_sequence_reader (const std::string& filename);


#endif
