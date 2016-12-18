#ifndef VITERBI_CPP
#define VITERBI_CPP

#include <string>
#include <vector>
#include "viterbi.h"

using namespace std;

//-------------------------------- constructors --------------------------------
/*Viterbi::Viterbi():Viterbi(1,1,"")
{
  //k = 1;
  //nstates = 1;
  //nsteps = 1;
  //viterbi = vector<vector<Probability> >(nstates,vector<Probability>(nsteps,Probability(1.0)));
}

Viterbi::Viterbi(int states)//:k(1),nstates(states)
{
  k = 1;
  nstates = states;
  nsteps = 1;
  viterbi = vector<vector<Probability> >(nstates,vector<Probability>(nsteps,Probability(1.0)));
}

Viterbi::Viterbi(int states, int kmer)
{
  k = kmer;
  nstates = states;
  nsteps = 1;
  viterbi = vector<vector<Probability> >(nstates,vector<Probability>(nsteps,Probability(1.0)));
}

Viterbi::Viterbi(int states, string obs)
{
  k = 1;
  nstates = states;
  nsteps = obs.length();
  viterbi = vector<vector<Probability> >(nstates,vector<Probability>(nsteps,Probability(1.0/nsteps)));
}
*/
Viterbi::Viterbi(int states, int kmer, string obs):k(kmer),nstates(states)
{
  //k = kmer;
  //nstates = states;
  nsteps = obs.empty()?0:obs.length()-k+1;
  viterbi = vector<vector<Probability> >(states,vector<Probability>(nsteps,Probability(0.0)));
}
//--------------------------------- operators ----------------------------------
vector<Probability>& Viterbi::operator[](int state)
{
  return viterbi[state];
}
//--------------------------------- functions ----------------------------------
void Viterbi::print_vertical()
{
  cout<<"viterbi: "<<endl;
  for(int i=0; i<nstates; i++)
  {
    cout<<'\t'<<i;
  }
  cout<<endl;
  for(int i=0; i<nsteps; i++)
  {
    cout<<i;
    for(int j=0; j<nstates; j++)
    {
      cout<<'\t'<<viterbi[j][i];
    }
    cout<<endl;
  }
}

void Viterbi::print_horizontal()
{
  cout<<"viterbi: "<<endl;
  for(int i=0; i<nsteps; i++)
  {
    cout<<'\t'<<i;
  }
  cout<<endl;
  for(int i=0; i<nstates; i++)
  {
    cout<<i;
    for(int j=0; j<nsteps; j++)
    {
      cout<<'\t'<<viterbi[i][j];
    }
    cout<<endl;
  }
}

#endif //VITERBI_CPP
