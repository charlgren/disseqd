#ifndef EMISSION_H
#define EMISSION_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include "probability.h"

using namespace std;

class Emission
{
public:
//-------------------------------- constructors --------------------------------
  /*Emission(Probability init = Probability(1.0));
  Emission(int states, Probability init = Probability(1.0));
  Emission(int states, int k, Probability init = Probability(1.0));*/
  Emission(int states=1, int k=1, Probability init = Probability(1.0));
//--------------------------------- operators ----------------------------------
  map<string,Probability>& operator[](int state);
//--------------------------------- functions ----------------------------------
  void normalize();
  void normalize(int state);
  vector<string> get_emissions();
  vector<vector<Probability> > get_probabilities();
  vector<Probability> get_probabilities(int state);
  Probability get_min(){return p_min;}
  void print_vertical();
  void print_horizontal();
  int get_nemissions(){return nemissions;}
  void add(int state, Probability p);
private:
//--------------------------------- variables ----------------------------------
  vector<map<string,Probability> > emit_p;
  int nstates;
  int nemissions;
  Probability init_p;
  Probability p_min = Probability(1.0);
  string alphabet;
//--------------------------------- functions ----------------------------------
  void initialize_emissions(string alphabet, int k, Probability init);
  int initialize_emissions(string alphabet, string prefix, int n, int k, int idx, Probability init);
};

#endif //EMISSION_H
