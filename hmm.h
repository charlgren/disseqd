#ifndef HMM_H
#define HMM_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath> // INFINITY
#include "functions.cpp"
#include "probability.h"
#include "start.cpp"
#include "emission.cpp"
#include "transition.cpp"
#include "viterbi.cpp"

using namespace std;

//flags set by '--debug' and '--verbose'
extern int debug_flag;
extern int verbose_flag;

class Hmm
{
public:
//-------------------------------- constructors --------------------------------
  /*Hmm();
  Hmm(int states);
  Hmm(int states, int kmer);*/
  Hmm(int states=1, int kmer=1);
//--------------------------------- operators ----------------------------------

//--------------------------------- functions ----------------------------------
  string read_fasta(string infile);
  void read_emissions_fasta(int state, string infile);
  void print_kmer(){cout<<">kmer\n"<<k<<"\n";}
  void print_states(){cout<<">states\n"<<nstates<<"\n";}
  void print_start(){start_p.print();}
  void print_transitions(){trans_p.print();}
  void print_emissions(){emit_p.print_vertical();}
  string decode(string &obs);
  void train(string &obs);
  string get_viterbi_path() const {return viterbi_path;}
  void set_viterbi_path(string &p) {viterbi_path = p;}
  double get_viterbi_prob() const {return viterbi_prob;}
  void set_viterbi_prob(double p) {viterbi_prob = p;}
  void write_decoding(string file);
  void write_model(string file);
  void read_model(string file);
private:
//--------------------------------- variables ----------------------------------
  // int verbose_flag = 0;
  double viterbi_prob = -INFINITY;
  string viterbi_path;
  int T;
  int k;
  int nstates;
  Start start_p;
  Emission emit_p;
  Transition trans_p;
  Viterbi viterbi_p;
  vector<Probability> scaling_p;
  vector<Probability> scaling_b;
  vector<Probability> scaling_g;
  vector<vector<Probability> > forward_p;
  vector<vector<Probability> > backward_p;
  vector<vector<Probability> > gamma_p;
  vector< vector< vector<Probability> > > xi_p;
//--------------------------------- functions ----------------------------------
  void forward(string &obs);
  void backward(string &obs);
  void gamma(string &obs);
  void xi(string &obs);
  void update(string &obs);
  void update_start_p();
  void update_trans_p();
  void update_emit_p(string &obs);
};

#endif //HMM_H
