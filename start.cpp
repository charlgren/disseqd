#ifndef START_CPP
#define START_CPP

#include <iostream>
//#include <numeric>//accumulate
#include "start.h"

using namespace std;

//-------------------------------- constructors --------------------------------
/*Start::Start()
{
  nstates = 1;
  start_p = vector<Probability>(nstates,Probability(1.0));
}*/

Start::Start(int states)
{
  nstates = states;
  //start_p = vector<Probability>(states,Probability(1.0/states));
  double diag[] = {0.9,0.1};//{0.3,0.7};
  start_p = vector<Probability>(states);
  for(int i=0; i<states; i++)
  {
    start_p[i] = Probability(diag[i]);
  }
}

Start::Start(vector<Probability> start_p)
{
  start_p = start_p;
  nstates = start_p.size();
}
//--------------------------------- operators ----------------------------------
Probability& Start::operator[](int state)
{
  return start_p[state];
}
//--------------------------------- functions ----------------------------------
/*void Start::set_value(int state, double prob)
{
  start_p[state] = prob;
}*/

void Start::normalize()
{
  Probability sum = Probability(0.0);
  //cout<<"normalizing start: ";
  for(auto&& prob:start_p)
  {
    //cout<<'\t'<<prob.get_value();
    sum += prob;
  }
  //cout<<endl;
  //cout<<"scaled with: "<<sum.get_value()<<endl;
  for(auto&& prob:start_p)
  {
    prob/=sum;
  }
}

void Start::print()
{
  cout<<">start"<<endl;
  for(int i=0; i<nstates; i++)
  {
    cout<<i<<'\t'<<start_p[i]<<endl;
  }
}

#endif //START_CPP
