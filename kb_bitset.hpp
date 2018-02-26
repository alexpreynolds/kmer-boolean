#ifndef KMER_BOOLEAN_BITSET_H_
#define KMER_BOOLEAN_BITSET_H_

#include <vector>
#include <exception>
#include <stdexcept>
#include <map>

namespace kmer_boolean
{
  class KB_Bitset
  {
  private:
    std::vector<unsigned char> _bitset;
    int _k;
    int _nbytes;
    int _max_bits;

    inline int
    _byte_offset_from_index(const int& idx)
    {
      int o = 0;
      int m = idx;
      while (m > 0) {
	m -= 8;
	if (m >= 0) {
	  o++;
	}
      }
      return o;
    }

  public:
    typedef enum
      {
	MerFilterUndefined = 0,
	MerFilterPresent,
	MerFilterAbsent,
	MerFilterAll,
      } MerFilterType;
    
    std::vector<unsigned char>& bitset() { return _bitset; }
    size_t nbytes() const { return _nbytes; }

    std::map<unsigned char, int> create_fmap(void) {
      std::map<unsigned char, int> m;
      m['A'] = 0;
      m['C'] = 1;
      m['T'] = 2;
      m['G'] = 3;
      return m;
    } 
    std::map<unsigned char, int> fmap = create_fmap();

    std::map<int, unsigned char> create_rmap(void) {
      std::map<int, unsigned char> m;
      m[0] = 'A';
      m[1] = 'C';
      m[2] = 'T';
      m[3] = 'G';
      return m;
    }
    std::map<int, unsigned char> rmap = create_rmap();

    void 
    reserve_for_k(const int& k)
    {
      //
      // we need 4^k entries (bits) to store kmer presence/absence 
      // state, which therefore takes 4^k/8 bytes (unsigned chars)
      // or 2^(2k-3) bytes
      //
      try {
	if (k < 1) {
	  throw std::domain_error("Error: k must be positive, non-zero integer");
	}
	_k = k;
	_nbytes = 1 << (2*_k - 3);
	_max_bits = (k > 1) ? 8 : 4;
	_nbytes = (_nbytes > 1) ? _nbytes : 1;
	_bitset.reserve(_nbytes);
      } catch (const std::exception& e) {
	std::cout << e.what() << std::endl;
	std::terminate();
      }
    }

    inline void 
    set_all(const bool b) 
    {
      for (int bidx = 0; bidx < _nbytes; ++bidx) {
	_bitset[bidx] = (b) ? 0xff : 0x00; 
      } 
    }

    void
    set(const int& idx, const bool& v)
    {
      int x = v ? 1 : 0;
      int r = idx % _max_bits;
      int offset = _byte_offset_from_index(idx);
      _bitset[offset] ^= (-x ^ _bitset[offset]) & (1UL << r);
    }

    bool
    get(const int& idx)
    {
      int r = idx % _max_bits;
      int offset = _byte_offset_from_index(idx);
      return (_bitset[offset] >> r) & 1U;
    }

    void
    get_all(const MerFilterType ft)
    {
      std::vector<size_t> mer_inc(_k, 0);

      for (int byte = 0; byte < _nbytes; ++byte) {
	unsigned char bits = _bitset[byte];
#ifdef DEBUG_FLAG
	std::cout << "byte [" << byte << "]  bytes [" << byte_to_binary(bits) << "]" << std::endl;
#endif
	for (int bidx = 0; bidx < _max_bits; ++bidx) {
	  bool bit = (bits >> bidx) & 1U;
	  // print mer
	  if (((ft == MerFilterPresent) && bit) || ((ft == MerFilterAbsent) && !bit) || (ft == MerFilterAll)) {
	    for (std::vector<size_t>::const_reverse_iterator inc = mer_inc.rbegin(); inc != mer_inc.rend(); ++inc) {
	      std::cout << rmap[*inc];
	    }
	    if (((ft == MerFilterPresent) || (ft == MerFilterAll)) && bit) {
	      std::cout << " found" << std::endl;
	    }
	    else if (((ft == MerFilterAbsent) || (ft == MerFilterAll)) && !bit) {
	      std::cout << " not found" << std::endl;
	    }
	  }
	  // increment the lowest digit
	  mer_inc[0]++;
	  for (int ridx = 0; ridx < _k - 1; ++ridx) {
	    if (mer_inc[ridx] == rmap.size()) {
	      // reset digit and increment next digit
	      mer_inc[ridx] = 0;
	      mer_inc[ridx + 1]++;
	    }
	  }
	}
      }
    }

    const char *
    byte_to_binary(const int& byte)
    {
      static char binary[9];
      binary[8] = '\0';
      int z;
      for (z = 128; z > 0; z >>= 1) {
	std::strcat(binary, ((byte & z) == z) ? "1" : "0");
      }
      return binary;
    }

    KB_Bitset();
    ~KB_Bitset();
  };

  KB_Bitset::KB_Bitset() {
  }

  KB_Bitset::~KB_Bitset() {
  }
}

#endif // KMER_BOOLEAN_BITSET_H_
