#ifndef VITERBI_H
#define VITERBI_H

#include <string>
#include <vector>
#include "probability.h"

using namespace std;

class Viterbi
{
public:
//-------------------------------- constructors --------------------------------
  /*Viterbi();
  Viterbi(int states);
  Viterbi(int states, int kmer);
  Viterbi(int states, string obs);*/
  Viterbi(int states=1, int kmer=1, string obs="");
//--------------------------------- operators ----------------------------------
  vector<Probability>& operator[](int state);
//--------------------------------- functions ----------------------------------
  void print_horizontal();
  void print_vertical();

private:
//--------------------------------- variables ----------------------------------
  vector<vector<Probability> > viterbi;
  int k;
  int nstates;
  int nsteps;
  //shoud these reside here instead of in the hmm class?
  //string viterbi_path;
  //double viterbi_prob;
};

#endif //VITERBI_H
