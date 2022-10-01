#include <fstream>
#include <iostream>
#include <ctime>
#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

// CLASS + METHOD DECLARATIONS ****************************************************************************
class BloomFilter{
  private:
    boost::dynamic_bitset<> filter;
    uint64_t* hashes;
    uint64_t bloom_size;
    int hash_number;
  
  public:
    BloomFilter(uint64_t n, int nf);
    void add_value(uint64_t code);
    bool is_present(uint64_t code);
    void build_from_file(string file, int kmer_size);
    void printTab();
};



// FUNCTION DECLARATIONS **********************************************************************************
uint64_t xorshift64(uint64_t x);
void multihash(uint64_t x, uint64_t * hashes, uint64_t nb_hashes, uint64_t max_val);

char complement(char nucl);
string end_rev_compl(string started_rev_compl, string kmer, int l, int index);
string try_rev_compl(string kmer);

uint64_t str_to_code(string kmer);
uint64_t next_code(string &old_kmer, char new_char);

char next_char(ifstream &stream);

string random_word(int kmer_size);
void compute_random_requests(BloomFilter bf, int n_rand_tests, int kmer_size);
void compute_specific_requests(BloomFilter bf, string kmer);



// ****** MAIN ********************************************************************************************
int main(int argc, char* argv[]){
  // PARAMETERS / ARGS PARSING
  string file;
  int kmer_size;
  uint64_t filter_size;
  int n_hash_functions;
  int n_rand_tests;

  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produces help message")
      ("file", po::value(&file), "fasta file to map into the bloom filter")
      ("k", po::value(&kmer_size), "size of kmers <= 31")
      ("n", po::value(&filter_size), "size of bloom filter (in bits) <= 2^34")
      ("nf", po::value(&n_hash_functions), "number of applied hash functions <= 64")  
      ("r", po::value(&n_rand_tests), "number of random is_present tests")  
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help") || !vm.count ("file") || !vm.count ("k") || !vm.count ("n") || !vm.count ("nf") || !vm.count ("r")) {
      cerr << desc << "\n";
  return 1;
  } 

  //BLOOM FILTER BUILDING
  cout << "Building filter ..." << endl;
  BloomFilter bf(filter_size, n_hash_functions);
  bf.build_from_file(file, kmer_size);
  cout << "Filter ready" << endl;

  //BLOOM FILTER REQUESTS
  compute_random_requests(bf, n_rand_tests, kmer_size);
  compute_specific_requests(bf, "AAACTTACTT");
  
  return 0;
}



// CLASS METHOD DEFINITIONS *******************************************************************************

BloomFilter::BloomFilter(uint64_t n, int nf){
  //Builds the filter, inits the values in class
  filter.resize(n);
  bloom_size = n;
  hashes = new uint64_t[nf];
  hash_number = nf;
}

void BloomFilter::add_value(uint64_t code){
  /*Hashes a kmer code multiple times and set to 1 the 
  filter bits at the hashes values locations */
  multihash(code, hashes, hash_number, bloom_size-1);
  for (int i=0; i<hash_number; i++){
    filter.set(hashes[i]);
  }
}

bool BloomFilter::is_present(uint64_t code){
  /*Hashes a kmer code multiple times and check all the 
  filter bits at the hashes values locations (logical AND)*/
  multihash(code, hashes, hash_number, bloom_size-1);
  for (int i=0; i<hash_number; i++){
    if (! filter.test(hashes[i])){
      return false;
    }
  }
  return true;
}

void BloomFilter::build_from_file(string file, int kmer_size){
  /* open file and reads it char by char while reading each kmer
  to hash them and add them in the bloom filter */
  ifstream infile;
  infile.open(file, ios::in); //open file read mode
  infile.ignore (numeric_limits<streamsize>::max(), '\n' ); //skip first line in FASTA

  //first kmer
  string current_kmer = "";
  for (int i=0; i<kmer_size; i++){ //need k < length(file)
    current_kmer += next_char(infile);
  }
  add_value(str_to_code(current_kmer));

  char ch = next_char(infile);

  //all file content
  while (ch != EOF){
    add_value(next_code(current_kmer, ch)); //modifies current_kmer in-place
    ch = next_char(infile);
  }

  infile.close();
}

void BloomFilter::printTab(){
  //For debug
  //cout << filter << endl;
  cout << filter.count() << endl;
}



// FUNCTION DEFINITIONS ***********************************************************************************

/** Hash function that uses xor properties
 * @param x a 64-bits value to hash
 * @return hashed value
 */
uint64_t xorshift64(uint64_t x)
{
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  return x;
}


/** Generate multiple hashes of the same value using sequences of xorshifts.
 * @param x The value to hash
 * @param hashes An array already allocated to fill with hash values
 * @param nb_hashes The number of hash values needed
 * @param max_val The maximum number that a hashed value can be.
 */
void multihash(uint64_t x, uint64_t * hashes, uint64_t nb_hashes, uint64_t max_val) {
	// Init 64 bits hashes
  hashes[0] = xorshift64(x);
  for (uint64_t i=1 ; i<nb_hashes ; i++)
    hashes[i] = xorshift64(hashes[i-1]);

  for (uint64_t i=0 ; i<nb_hashes ; i++)
  	hashes[i] %= max_val + 1;
}


char complement(char nucl){
  //Returns the complement of a nucleotide
  switch (nucl){
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
  }
  return ' ';
}

string end_rev_compl(string started_rev_compl, string kmer, int l, int index){
  /* Used in try_rev_compl(), ends up the process of finding the 
  reverse complement of a DNA seq started in try_rev_compl() */
  for (int i=index; i >= 0; i--){
    started_rev_compl += complement(kmer[i]);
  }
  return started_rev_compl;
}

string try_rev_compl(string kmer){
  /* Checks the lexicographical order of a kmer and 
  its reverse complement then returns the first */
  char alt; 
  string rev_compl = "";
  int len = kmer.length();
  
  for (int i=0; i<len; i++){
    alt = complement(kmer[len-i-1]);

    if (kmer[i] > alt){ //case revcomp better
      rev_compl += alt;
      return end_rev_compl(rev_compl, kmer, len, len-i-2);
    }
    else if (kmer[i] < alt){ //case kmer better
      return kmer;
    }
    else { //case same letters, keep building rev_compl just in case
      rev_compl += alt; 
    }
  }
  return kmer;
}


uint64_t str_to_code(string kmer){ 
  // returns a hash of a kmer or of its rev complement
  kmer = try_rev_compl(kmer);
  uint64_t res = 0;

  for (int i = 0; i < kmer.size(); ++i) {
    switch (kmer[i]){
      case 'A':
        break;
      case 'C':
        res += 1 * pow(4, i);
        break;
      case 'G':
        res += 2 * pow(4, i);
        break;
      case 'T':
        res += 3 * pow(4, i);
        break;
    }
  }

  return res+1;
} 


uint64_t next_code(string &old_kmer, char new_char){
  // modifies old_kmer in-place and returns the hash of the modified kmer
  old_kmer.erase(0, 1);
  old_kmer += new_char;
  return str_to_code(old_kmer);
}


char next_char(ifstream &stream){
  //reads next character in file. Omits 'N'.
  char res;

  do {
    res = stream.get();
  }
  while (res == '\n' || res == 'N');

  return res;
}


void compute_random_requests(BloomFilter bf, int n_rand_tests, int kmer_size){
  //Generates random words and tests their absence in the bloom filter
  srand((unsigned int)time(NULL));
  uint64_t tmp = rand();

  static const char nucls[4] = {'A', 'T', 'C', 'G'};
  string random_word;
  
  int found_count = 0;

  for (int i=0; i<n_rand_tests; i++){
    random_word = "";
    for (int j=0; j<kmer_size; j++){
      tmp = xorshift64(tmp);
      random_word += nucls[tmp % 4];
    }

    if (bf.is_present(str_to_code(random_word))){
      found_count ++;
    }
  }

  cout << found_count << " random words founds over " << n_rand_tests << " tested." << endl;
}

void compute_specific_requests(BloomFilter bf, string kmer){
  //Test the absence of a kmer in the bloom filter
  string presence = bf.is_present(str_to_code(kmer)) ? " has been found." : " is absent.";
  cout << "\"" << kmer << "\"" << presence << endl;
}