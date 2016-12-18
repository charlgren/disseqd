#ifndef START_H
#define START_H

#include <iostream>
#include "probability.h"

using namespace std;

class Start
{
public:
//-------------------------------- constructors --------------------------------
  //Start();
  Start(int states=1);
  Start(vector<Probability> start_p);
//--------------------------------- operators ----------------------------------
  Probability& operator[](int state);
//--------------------------------- functions ----------------------------------
  //void set_value(int state, double prob);
  void print();
  void normalize();

private:
//--------------------------------- variables ----------------------------------
  vector<Probability> start_p;
  int nstates;
};

#endif //START_H
