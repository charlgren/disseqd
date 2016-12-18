#ifndef EMISSION_CPP
#define EMISSION_CPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include "emission.h"

using namespace std;

//-------------------------------- constructors --------------------------------
/*Emission::Emission(Probability init)
{
  nstates = 1;
  alphabet = "ACGTN";
  nemissions = 1;
  emit_p = vector<map<string,Probability> >(nstates,map<string,Probability>());
  init_p = init;
  initialize_emissions(alphabet,1,init);
}

Emission::Emission(int states, Probability init)
{
  nstates = states;
  alphabet = "ACGTN";
  nemissions = pow(alphabet.length(),1);
  emit_p = vector<map<string,Probability> >(states,map<string,Probability>());
  init_p = init;
  initialize_emissions(alphabet,1,init);
}

Emission::Emission(int states, int k, Probability init)
{
  nstates = states;
  alphabet = "ACGTN";
  nemissions = pow(alphabet.length(),k);
  emit_p = vector<map<string,Probability> >(states,map<string,Probability>());
  init_p = init;
  initialize_emissions(alphabet,k,init);
}*/

Emission::Emission(int states, int k, Probability init):nstates(states)
{
  alphabet = "ACGTN";
  nemissions = pow(alphabet.length(),k);
  emit_p = vector<map<string,Probability> >(states,map<string,Probability>());
  init_p = init;
  initialize_emissions(alphabet,k,init);
}
//--------------------------------- operators ----------------------------------
map<string,Probability>& Emission::operator[](int state)
{
  return emit_p[state];
}
//--------------------------------- functions ----------------------------------
void Emission::add(int state, Probability p)
{
  //for(int i=0; i<nstates; i++)
  //{
    for(auto&& pair:emit_p[state])
    {
      pair.second+=p;
    }
  //}
}

void Emission::initialize_emissions(string alphabet, int k, Probability init)
{
  int idx = initialize_emissions(alphabet,"",alphabet.length(),k,0,init);
}

//don't need idx without counts
int Emission::initialize_emissions(string alphabet, string prefix, int n, int k, int idx, Probability init)
{
  if(k==0)
  {
    for(int s=0; s<nstates; s++)
    {
      emit_p[s][prefix]=init;
      //emit_p[s][prefix]=Probability(0.0000000001);
      //emit_p[s][prefix]=Probability(1.0/nemissions);
      //emit_p[s][prefix]=Probability(1.0/alphabet.length());
      //emit_p[s][prefix]=Probability(1.0);
    }
//    cout<<"adding: "<<prefix<<" at: "<<idx<<endl;
    return ++idx;
  }
  for(int i=0; i<n; ++i)
  {
    string newprefix = prefix + alphabet.at(i);
    idx = Emission::initialize_emissions(alphabet, newprefix, n, k-1, idx, init);
  }
  return idx;
}

void Emission::normalize()
{
  for(int state=0; state<nstates; state++)
  {
    Emission::normalize(state);
  }
}

void Emission::normalize(int state)
{
  vector<Probability> vsum = vector<Probability>(nemissions,Probability(0.0));
  //Probability sum = Probability(0.0);
  int i = 0;
  for(auto&& pair:emit_p[state])
  {
    //cerr<<pair.first<<i/alphabet.length()<<endl;
    vsum[i/alphabet.length()] += pair.second;
    //vsum[i/nemissions] += pair.second;
    i++;
    //sum += pair.second;
  }
  //Probability sum2 = Probability(0.0);
  i = 0;
  for(auto&& pair:emit_p[state])
  {
    //cerr<<state<<'\t'<<pair.first<<'\t'<<vsum[i/alphabet.length()]<<endl;
    pair.second/=vsum[i/alphabet.length()];
    //pair.second/=vsum[i/nemissions];
    i++;
    p_min = pair.second<p_min?pair.second:p_min;
    //pair.second/=sum;
    //sum2 += pair.second;
  }
  //cout<<"normalizing state: "<<state<<"\tfrom: "<<sum.get_value()<<"\tto: "<<sum2.get_value()<<endl;
}

vector<string> Emission::get_emissions()
{
  vector<string> emissions;
  for(auto&& pair:emit_p[0])
  {
    emissions.push_back(pair.first);
  }
  return emissions;
}

vector<vector<Probability> > Emission::get_probabilities()
{
  vector<string> emissions = get_emissions();
  vector<vector<Probability> > probabilities(nstates,vector<Probability>(nemissions));
  for(int i=0; i<nstates; i++)
  {
    probabilities[i] = get_probabilities(i);
  }
  return probabilities;
}

vector<Probability> Emission::get_probabilities(int state)
{
  vector<string> emissions = get_emissions();
  vector<Probability> probabilities(nemissions);
  for(int i=0; i<nemissions; i++)
  {
    probabilities[i] = emit_p[state][emissions[i]];
  }
  return probabilities;
}

void Emission::print_vertical()
{
  cout<<">emissions"<<endl;
  vector<string> emissions = get_emissions();
  vector<vector<Probability> > probabilities = get_probabilities();

  vector<Probability> sum = vector<Probability>(nstates,Probability(0.0));
  for(int i=0; i<nstates; i++)
  {
    cout<<'\t'<<i;
  }
  cout<<endl;
  for(int i=0; i<nemissions; i++)
  {
    for(int j=0; j<nstates; j++)
    {
      if(j==0)
      {
        cout<<emissions[i];
      }
      sum[j]  += probabilities[j][i];
      cout<<'\t'<<probabilities[j][i];
    }
    cout<<endl;
  }
/*  cout<<"sum:\t";
  for(auto&& state_sum:sum)
  {
    cout<<'\t'<<state_sum;
  }
  cout<<endl;
*/
}

void Emission::print_horizontal()
{
  vector<string> emissions = get_emissions();
  for(int i=0; i<nemissions; i++)
  {
    cout<<'\t'<<emissions[i];
  }
  cout<<endl;
  for(int i=0; i<nstates; i++)
  {
    cout<<i;
    vector<Probability> probabilities = get_probabilities(i);
    for(int j=0; j<nemissions; j++)
    {
      cout<<'\t'<<probabilities[j];
    }
    cout<<endl;
  }
}

#endif //EMISSION_CPP
