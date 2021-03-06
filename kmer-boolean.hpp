#ifndef KMER_BOOLEAN_H_
#define KMER_BOOLEAN_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif /* getline() support */

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <exception>
#include <getopt.h>
#include <iostream>
#include "kb_bitset.hpp"

#define KMER_BOOLEAN_LINE_MAX 268435456
#define KMER_BOOLEAN_STDIN_BUF_MAX 8192
#define KMER_BOOLEAN_SEQ_BUF_MAX KMER_BOOLEAN_STDIN_BUF_MAX*1024

namespace kmer_boolean
{
  class KB 
  {
  private:
    int _k;
    bool _read_in_all_sequences_at_once;
    std::string _query_kmer;
    KB_Bitset _bitset;
    KB_Bitset::MerFilterType _filter_type;
    
  public:
    std::vector<std::string> sequences;

    static const std::string client_name;
    static const std::string client_version;
    static const std::string client_authors;
    
    void process_sequences_by_chunks(void);
    void process_sequence_chunk_buffer(std::vector<char>& _seq_chunk_buf_vec, std::streamsize& _seq_chunk_nchars);
    void read_all_sequences(void);
    void write_sequences(void);
    void initialize_bitset(void);
    void process_all_sequences(void);
    void get_mers_with_state(void);
    void test_mer(void);

    void initialize_command_line_options(int argc, char** argv);
    void print_usage(FILE* os);
    void print_version(FILE* os);
    
    std::string client_kmer_boolean_opt_string(void);
    struct option* client_kmer_boolean_long_options(void);
    std::string client_kmer_boolean_name(void);
    std::string client_kmer_boolean_version(void);
    std::string client_kmer_boolean_authors(void);
    std::string client_kmer_boolean_usage(void);
    std::string client_kmer_boolean_description(void);
    std::string client_kmer_boolean_io_options(void);
    std::string client_kmer_boolean_general_options(void);

    int k() const { return _k; }
    void k(const int& k) { _k = k; }
    
    bool read_in_all_sequences_at_once() const { return _read_in_all_sequences_at_once; }
    void read_in_all_sequences_at_once(const bool& r) { _read_in_all_sequences_at_once = r; }

    std::string query_kmer() const { return _query_kmer; }
    void query_kmer(const std::string& qk) { _query_kmer = qk; }

    KB_Bitset& bitset() { return _bitset; }

    KB_Bitset::MerFilterType filter_type() const { return _filter_type; }
    void filter_type(const KB_Bitset::MerFilterType& mft) { _filter_type = mft; }

    static KB_Bitset::MerFilterType& default_filter_type() { 
      static KB_Bitset::MerFilterType dft = KB_Bitset::MerFilterPresent;
      return dft;
    }

    KB();
    ~KB();
  };

  KB::KB() {
    k(-1);
    query_kmer("");
    filter_type(KB_Bitset::MerFilterUndefined);
  }
  
  KB::~KB() {
  }
}

#endif // KMER_BOOLEAN_H_
