#include "kmer-boolean.hpp"

const std::string kmer_boolean::KB::client_name = "kmer-boolean";
const std::string kmer_boolean::KB::client_version = "1.1";
const std::string kmer_boolean::KB::client_authors = "Alex Reynolds";

int
main(int argc, char** argv)
{
  kmer_boolean::KB kb;

  // setup
  kb.initialize_command_line_options(argc, argv);
  kb.initialize_bitset();
  std::ios_base::sync_with_stdio(false);
  
  // parse sequence input
  if (kb.read_in_all_sequences_at_once()) {
    kb.read_all_sequences();
    kb.process_all_sequences();
  }
  else {
    try {
      kb.process_sequences_by_chunks();
    }
    catch (const std::invalid_argument& ia) {
      std::cerr << ia.what() << std::endl;
      return EXIT_FAILURE;
    } 
  }
  // write kmer report
  if (kb.query_kmer().empty()) {
    kb.get_mers_with_state();
  }
  else {
    kb.test_mer();
  }
  
  return EXIT_SUCCESS;
}

void
kmer_boolean::KB::process_sequences_by_chunks(void)
{
  if (this->k() > KMER_BOOLEAN_SEQ_BUF_MAX) {
    throw std::invalid_argument("Error: k is larger than intermediate sequence buffer size");
  }
  bool _within_fasta_record_header = false;
  std::vector<char> _stdin_buf_vec(KMER_BOOLEAN_STDIN_BUF_MAX, 0);
  std::vector<char> _seq_buf_vec(KMER_BOOLEAN_SEQ_BUF_MAX, 0);
  std::streamsize _stdin_nchars = 0;
  std::streamsize _seq_buf_size = 0;
  while (std::cin.read(_stdin_buf_vec.data(), _stdin_buf_vec.size())) {
    _stdin_nchars = std::cin.gcount();
#ifdef DEBUG_FLAG
    std::string _buf_str(_stdin_buf_vec.begin(), _stdin_buf_vec.end());
    std::cerr << "read in " << _stdin_nchars << " chars from stdin: [" << _buf_str << "]" << std::endl;
#endif
    /*
      If buffer starts or contains '>' character, this starts a FASTA 
      record header, and a newline subsequently starts a FASTA record 
      sequence.
    */
    for (std::streamsize i = 0; i < _stdin_nchars; ++i) {
      if (_stdin_buf_vec[i] == '>') {
        _within_fasta_record_header = true;
        /*
          If _seq_buf_size is greater than zero, we need to process 
          whatever is in _seq_buf_vec and reset the size parameter
        */
        this->process_sequence_chunk_buffer(_seq_buf_vec, _seq_buf_size);
        _seq_buf_size = 0;
      }
      else if ((_stdin_buf_vec[i] == '\n') && _within_fasta_record_header) {
        _within_fasta_record_header = false;
      }
      else if ((_stdin_buf_vec[i] != '\n') && (_stdin_buf_vec[i] != 'N') && !_within_fasta_record_header) {
#ifdef DEBUG_FLAG
        std::fprintf(stderr, "%c | %c\n", _stdin_buf_vec[i], (_within_fasta_record_header ? 'T' : 'F'));
#endif
        _seq_buf_vec[_seq_buf_size++] = std::toupper(_stdin_buf_vec[i]);
        if (_seq_buf_size == KMER_BOOLEAN_SEQ_BUF_MAX) {
          this->process_sequence_chunk_buffer(_seq_buf_vec, _seq_buf_size);
          // copy k-1 chars from end of _seq_buf_vec and reset _seq_buf_size
          for (std::streamsize j = 0; j < this->k() - 1; ++j) {
            std::streamsize offset = KMER_BOOLEAN_SEQ_BUF_MAX - this->k() + j + 1;
            _seq_buf_vec[j] = _seq_buf_vec[offset];
#ifdef DEBUG_FLAG
            std::fprintf(stderr, " (writing [%c] from offset [%zu] to offset [%zu]\n", _seq_buf_vec[j], offset, j);
#endif
          }
          _seq_buf_size = this->k() - 1;
        }
#ifdef DEBUG_FLAG
        else {
          std::fprintf(stderr, "current state ->\n");
          for (std::streamsize j = 0; j < _seq_buf_size; ++j) {
            std::fprintf(stderr, " [%c]\n", _seq_buf_vec[j]);
          }
        }
#endif
      }
    }
  }
  _stdin_nchars = std::cin.gcount();
#ifdef DEBUG_FLAG
  std::string _end_buf_str(_stdin_buf_vec.begin(), _stdin_buf_vec.end());
  std::cerr << "read in " << _stdin_nchars << " chars from stdin: [" << _end_buf_str << "]" << std::endl;
#endif
  for (std::streamsize i = 0; i < _stdin_nchars; ++i) {
    if (_stdin_buf_vec[i] == '>') {
      _within_fasta_record_header = true;
      /*
        If _seq_buf_size is greater than zero, we need to process 
        whatever is in _seq_buf_vec and reset the size parameter
      */
      this->process_sequence_chunk_buffer(_seq_buf_vec, _seq_buf_size);
      _seq_buf_size = 0;
    }
    else if ((_stdin_buf_vec[i] == '\n') && _within_fasta_record_header) {
      _within_fasta_record_header = false;
    }
    else if ((_stdin_buf_vec[i] != '\n') && (_stdin_buf_vec[i] != 'N') && !_within_fasta_record_header) {
#ifdef DEBUG_FLAG
      std::fprintf(stderr, "%c | %c\n", _stdin_buf_vec[i], (_within_fasta_record_header ? 'T' : 'F'));
#endif
      _seq_buf_vec[_seq_buf_size++] = std::toupper(_stdin_buf_vec[i]);
    }
  }
#ifdef DEBUG_FLAG
  std::fprintf(stderr, "end _seq_buf_size [%zu]\n", _seq_buf_size);
#endif
  this->process_sequence_chunk_buffer(_seq_buf_vec, _seq_buf_size);
}

void
kmer_boolean::KB::process_sequence_chunk_buffer(std::vector<char>& _seq_chunk_buf_vec, std::streamsize& _seq_chunk_nchars) {
  // if the sequence chunk buffer is shorter than k, then we do not process it
  if (_seq_chunk_nchars < this->k()) {
    return;
  }
#ifdef DEBUG_FLAG
  std::string _buf_str;
  std::streamsize _pos = 0;
  for (std::vector<char>::const_iterator i = _seq_chunk_buf_vec.begin(); _pos < _seq_chunk_nchars; ++i) {
    _buf_str += *i;
    ++_pos;
  }
  std::fprintf(stderr, " --> scanning [%s] for %dmers\n", _buf_str.c_str(), this->k());
#endif
  std::deque<char> window(_seq_chunk_buf_vec.begin(), _seq_chunk_buf_vec.begin() + this->k());
  for (std::streamsize i = this->k(); i <= _seq_chunk_nchars; ++i) {
    std::string mer(window.begin(), window.end());
#ifdef DEBUG_FLAG
    std::fprintf(stderr, " ----> bitset [%s]\n", mer.c_str());
#endif
    int idx = 0;
    for (int j = this->k() - 1; j > -1; --j) {
      idx += bitset().fmap[mer[j]] * 1 << (2*(this->k() - 1 - j));
    }
    bitset().set(idx, true);
    window.pop_front();
    window.push_back(_seq_chunk_buf_vec[i]);
  }
}

void
kmer_boolean::KB::get_mers_with_state(void)
{
  try {
    switch (this->filter_type()) {
    case (KB_Bitset::MerFilterPresent):
    case (KB_Bitset::MerFilterAbsent):
    case (KB_Bitset::MerFilterAll):
      this->bitset().get_all(this->filter_type());
      break;
    case (KB_Bitset::MerFilterUndefined):
    default:
      throw std::invalid_argument("Error: mer filter type is undefined");
    }
  } catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    std::terminate();
  }
}

void
kmer_boolean::KB::test_mer(void)
{
  int idx = 0;
  for (int j = this->k() - 1; j > -1; --j) {
    idx += bitset().fmap[this->query_kmer()[j]] * 1 << (2*(this->k() - 1 - j));
  }
  bool test_result = bitset().get(idx);
  std::cout << this->query_kmer() << ((test_result) ? " found" : " not found") << std::endl;
}

void
kmer_boolean::KB::process_all_sequences(void)
{
  for (std::vector<std::string>::const_iterator seq = this->sequences.begin(); seq != this->sequences.end(); ++seq) {
    if ((*seq).length() < (unsigned long)this->k()) {
      std::ostringstream err;
      err << "Error: k (" << this->k() << ") is longer than the input sequence length (" << (*seq).length() << ")";
      throw std::invalid_argument(err.str());
    }
    std::deque<unsigned char> window((*seq).begin(), (*seq).begin() + this->k());
    for (size_t i = this->k(); i <= (*seq).length(); ++i) {
      std::string mer(window.begin(), window.end());
      window.pop_front();
      window.push_back((*seq)[i]);
      int idx = 0;
      for (int j = this->k() - 1; j > -1; --j) {
        idx += bitset().fmap[mer[j]] * 1 << (2*(this->k() - 1 - j));
      }
      bitset().set(idx, true);
    }
  }
}

void
kmer_boolean::KB::initialize_bitset(void)
{
  bitset().reserve_for_k(k());
  bitset().set_all(false);
}

void
kmer_boolean::KB::read_all_sequences(void)
{
  std::string sequence("");
  for (std::string line; std::getline(std::cin, line);) {
    if (line.find(">") == 0) {
      if (!sequence.empty()) {
        this->sequences.push_back(sequence);
        sequence = "";
      }
    }
    else {
      line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
      sequence += line;
    }
  }
  this->sequences.push_back(sequence);
}

void
kmer_boolean::KB::write_sequences(void)
{
  for (std::vector<std::string>::const_iterator seq = this->sequences.begin(); seq != this->sequences.end(); ++seq) {
    std::cout << *seq << std::endl; 
  }
}

std::string
kmer_boolean::KB::client_kmer_boolean_opt_string(void)
{
  static std::string _s("k:q:palrhv?");
  return _s;
}

struct option*
kmer_boolean::KB::client_kmer_boolean_long_options(void)
{
  static struct option _k = { "k",                   required_argument,   NULL,    'k' };
  static struct option _q = { "query-kmer",          required_argument,   NULL,    'q' };
  static struct option _p = { "present",             no_argument,         NULL,    'p' };
  static struct option _a = { "absent",              no_argument,         NULL,    'a' };
  static struct option _l = { "all",                 no_argument,         NULL,    'l' };
  static struct option _r = { "read-in-all-at-once", no_argument,         NULL,    'r' };
  static struct option _h = { "help",                no_argument,         NULL,    'h' };
  static struct option _v = { "version",             no_argument,         NULL,    'v' };
  static struct option _0 = { NULL,                  no_argument,         NULL,     0  };
  static std::vector<struct option> _options;
  _options.push_back(_k);
  _options.push_back(_q);
  _options.push_back(_p);
  _options.push_back(_a);
  _options.push_back(_l);
  _options.push_back(_r);
  _options.push_back(_h);
  _options.push_back(_v);
  _options.push_back(_0);
  return &_options[0];
}

void
kmer_boolean::KB::initialize_command_line_options(int argc, char** argv)
{
  int client_long_index;
  int client_opt = getopt_long(argc,
                               argv,
                               this->client_kmer_boolean_opt_string().c_str(),
                               this->client_kmer_boolean_long_options(),
                               &client_long_index);
  int _k = -1;
  this->read_in_all_sequences_at_once(false);

  opterr = 0; /* disable error reporting by GNU getopt */

  bool filter_set = false;

  while (client_opt != -1) {
    switch (client_opt) {
    case 'k':
      std::sscanf(optarg, "%d", &_k);
      this->k(_k);
      break;
    case 'q':
      this->query_kmer(optarg);
      break;
    case 'p':
      this->filter_type(KB_Bitset::MerFilterPresent);
      filter_set = true;
      break;
    case 'a':
      this->filter_type(KB_Bitset::MerFilterAbsent);
      filter_set = true;
      break;
    case 'l':
      this->filter_type(KB_Bitset::MerFilterAll);
      filter_set = true;
      break;
    case 'r':
      this->read_in_all_sequences_at_once(true);
      break;
    case 'h':
      this->print_usage(stdout);
      std::exit(EXIT_SUCCESS);
    case 'v':
      this->print_version(stdout);
      std::exit(EXIT_SUCCESS);
    case '?':
      this->print_usage(stdout);
      std::exit(EXIT_SUCCESS);
    default:
      break;
    }
    client_opt = getopt_long(argc,
                             argv,
                             this->client_kmer_boolean_opt_string().c_str(),
                             this->client_kmer_boolean_long_options(),
                             &client_long_index);
  }

  bool error_flagged = false;
  
  if (this->k() == -1) {
    std::fprintf(stderr, "Error: Specify k value\n");
    error_flagged = true;
  }

  if (!this->query_kmer().empty() && this->query_kmer().length() != (size_t)this->k()) {
    std::fprintf(stderr, "Error: Specify query kmer with length equivalent to k\n");
    error_flagged = true;
  }

  if (error_flagged) {
    this->print_usage(stderr);
    std::exit(ENODATA);
  }

  if (this->query_kmer().empty() && !filter_set) {
    this->filter_type(kmer_boolean::KB::default_filter_type());
  }
}

std::string
kmer_boolean::KB::client_kmer_boolean_name(void)
{
  static std::string _s(kmer_boolean::KB::client_name);
  return _s;
}

std::string
kmer_boolean::KB::client_kmer_boolean_version(void)
{
  static std::string _s(kmer_boolean::KB::client_version);
  return _s;
}

std::string
kmer_boolean::KB::client_kmer_boolean_authors(void)
{
  static std::string _s(kmer_boolean::KB::client_authors);
  return _s;
}

std::string
kmer_boolean::KB::client_kmer_boolean_usage(void)
{
  static std::string _s("\n"						\
                        "  Usage:\n"		\
                        "\n"						\
                        "  $ kmer_boolean [arguments] < input\n");
  return _s;
}

std::string
kmer_boolean::KB::client_kmer_boolean_description(void)
{
  static std::string _s("  Test if specified kmer is or is not in a set of\n" \
                        "  FASTA sequences provided on standard input, for a\n" \
                        "  given k, returning the according 'true' or 'false'\n" \
                        "  result.\n");
  return _s;
}

std::string
kmer_boolean::KB::client_kmer_boolean_io_options(void)
{
  static std::string _s("  General Options:\n\n"              \
                        "  --k=n                          K-value for kmer length (integer)\n" \
                        "  --query-kmer=s                 Test or query kmer (string)\n" \
                        " [--present | --absent | --all]  If query kmer is omitted, print all\n" \
                        "                                 present or absent kmers, or all kmers\n" \
                        "  --read-in-all-at-once          Read all sequence data into memory\n" \
                        "                                 before processing (not recommended for\n" \
                        "                                 very large FASTA inputs; default is to\n" \
                        "                                 stream through input in chunks)\n");
  return _s;
}

std::string
kmer_boolean::KB::client_kmer_boolean_general_options(void)
{
  static std::string _s("  Process Flags:\n\n"				\
                        "  --help                         Show this usage message\n" \
                        "  --version                      Show binary version\n");
  return _s;
}

void
kmer_boolean::KB::print_usage(FILE* os)
{
  std::fprintf(os,
               "%s\n"						     \
               "  version: %s\n"				     \
               "  author:  %s\n"				     \
               "%s\n"						     \
               "%s\n"						     \
               "%s\n"						     \
               "%s\n",
               this->client_kmer_boolean_name().c_str(),
               this->client_kmer_boolean_version().c_str(),
               this->client_kmer_boolean_authors().c_str(),
               this->client_kmer_boolean_usage().c_str(),
               this->client_kmer_boolean_description().c_str(),
               this->client_kmer_boolean_io_options().c_str(),
               this->client_kmer_boolean_general_options().c_str());
}

void
kmer_boolean::KB::print_version(FILE* os)
{
  std::fprintf(os,
               "%s\n"						     \
               "  version: %s\n"				     \
               "  author:  %s\n",
               this->client_kmer_boolean_name().c_str(),
               this->client_kmer_boolean_version().c_str(),
               this->client_kmer_boolean_authors().c_str());
}
