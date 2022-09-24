#include <stdio.h>
#include <stdint.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>

using namespace std;
using byte = unsigned char;

class BloomFilter{
  public:
    byte filter;

    BloomFilter(int n, int nf){
      byte filter[n/8 + 1] = {'0', '0', '0', '0', '0', '0', '0', '0'};
      cout << filter;
    }
};

char nextChar(ifstream &stream){
  char res;

  do {
    res = stream.get();
  }
  while (res == '\n' || res == 'N');

  return res;
}

uint64_t strToCode(string kmer){ 
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

uint64_t nextCode(uint64_t oldCode, char oldChar, char newChar, int k){
  switch (oldChar){
      case 'A':
        oldCode -= 1;
        break;
      case 'C':
        oldCode -= 2;
        break;
      case 'G':
        oldCode -= 3;
        break;
      case 'T':
        oldCode -= 4;
        break;
    }

    oldCode /= 4;

    switch (newChar){
      case 'A':
        break;
      case 'C':
        oldCode += 1 * pow(4, k-1);
        break;
      case 'G':
        oldCode += 2 * pow(4, k-1);
        break;
      case 'T':
        oldCode += 3 * pow(4, k-1);
        break;
    }
  return oldCode+1;
}


// ****** MAIN ****************************************************************
int main(){
  ifstream infile;
  infile.open("cov.fasta", ios::in); //open file read mode
  infile.ignore (numeric_limits<streamsize>::max(), '\n' ); //skip first line in FASTA

  char ch = nextChar(infile);
  char oldCh;
  int64_t code = strToCode("AA");

  while (ch != EOF){
    cout << ch << endl;
    oldCh = ch;
    ch = nextChar(infile);
    code = nextCode(code, oldCh, ch, 2);
    cout << code << endl;
  }

  cout << strToCode("AA") << endl;

  BloomFilter bf(32, 8);
  cout << bf.filter << endl;
  
  infile.close();
  return 0;
}