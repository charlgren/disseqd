#ifndef TRANSITION_H
#define TRANSITION_H

#include <string>
#include <vector>
#include <iostream>
#include "probability.h"

using namespace std;

class Transition
{
public:
//-------------------------------- constructors --------------------------------
  Transition(int states = 1);
//--------------------------------- operators ----------------------------------
  vector<Probability>& operator[](int state);
//--------------------------------- functions ----------------------------------
  void normalize();
  void normalize(int state);
//  void add_transition_probabilities(int from_state, valarray<double> probs){trans_p[from_state]=probs;}
//  double get_transition_probability(int from_state, int to_state){return trans_p[from_state][to_state];}
//  void set_transition_probability(int from_state, int to_state, double prob){trans_p[from_state][to_state]=prob;}
  void print();
private:
//--------------------------------- variables ----------------------------------
  vector<vector<Probability> > trans_p;
  int nstates;
};

#endif //TRANSITION_H
