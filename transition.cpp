#ifndef TRANSITION_CPP
#define TRANSITION_CPP

#include <string>
#include <vector>
#include <iostream>
#include "transition.h"

using namespace std;

//-------------------------------- constructors --------------------------------
Transition::Transition(int states):nstates(states)
{
  trans_p = vector<vector<Probability> >(nstates,vector<Probability>(nstates,Probability(1.0/nstates)));
  double diag[] = {1.-1./50,1.-1./5};//{0.99,0.21};//{0.995,0.55};//0.99(5),0.55
  vector<double> diag_p(diag, diag + sizeof(diag) / sizeof(double));
  for(int i=0; i<nstates; i++)
  {
    for(int j=0; j<nstates; j++)
    {
      trans_p[i][j]=(i==j?Probability(diag_p[i]):Probability((1.0-diag_p[i])/(nstates-1)));
    }
  }
}
//--------------------------------- operators ----------------------------------
vector<Probability>& Transition::operator[](int state)
{
  return trans_p[state];
}
//--------------------------------- functions ----------------------------------
void Transition::normalize()
{
  for(int state=0; state<nstates; state++)
  {
    normalize(state);
  }
}

void Transition::normalize(int state)
{
  Probability sum = Probability(0.0);
  //cout<<"normalizing trans: ";
  for(auto&& prob:trans_p[state])
  {
    //cout<<'\t'<<prob.get_value();
    sum += prob;
  }
  //cout<<endl;
  Probability sum2 = Probability(0.0);
  for(auto&& prob:trans_p[state])
  {
    prob/=sum;
    sum2 += prob;
  }
  //cout<<"normalizing state: "<<state<<"\tfrom: "<<sum.get_value()<<"\tto: "<<sum2.get_value()<<endl;
}

void Transition::print()
{
  cout<<">transitions"<<endl;
  for(int i=0; i<nstates; i++)
  {
    cout<<'\t'<<i;
  }
  cout<<endl;
  for(int i=0; i<nstates; i++)
  {
    cout<<i;
    for(int j=0; j<nstates; j++)
    {
      cout<<'\t'<<trans_p[i][j];
    }
    cout<<endl;
  }
}

#endif //TRANSITION_CPP
