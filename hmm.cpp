#ifndef HMM_CPP
#define HMM_CPP

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>
#include "hmm.h"

int debug_flag;
int verbose_flag;

using namespace std;

//-------------------------------- constructors --------------------------------
/*Hmm::Hmm()
{
  k = 1;
  nstates = 1;
  start_p = Start();
  emit_p = Emission();
  trans_p = Transition();
  viterbi_p = Viterbi();
}

Hmm::Hmm(int states)
{
  k = 1;
  nstates = states;
  start_p = Start(states);
  emit_p = Emission(states);
  trans_p = Transition(states);
  viterbi_p = Viterbi(states);
}

Hmm::Hmm(int states, int kmer)
{
  k = kmer;
  nstates = states;
  start_p = Start(states);
  emit_p = Emission(states,kmer);
  trans_p = Transition(states);
  viterbi_p = Viterbi(states,kmer);
}*/

Hmm::Hmm(int states, int kmer):k(kmer),nstates(states)
{
  start_p = Start(states);
  emit_p = Emission(states,kmer);
  trans_p = Transition(states);
  viterbi_p = Viterbi(states,kmer);
}
//--------------------------------- functions ----------------------------------
void Hmm::train(string &obs)
{
  T = obs.length()-k+1;//+1//size of forward, backward, gamma = length-k+1+1
  //cout<<"T="<<T<<endl;
  forward(obs);
  backward(obs);
  gamma(obs);
  xi(obs);
  update(obs);
}

void Hmm::forward(string &obs)
{//verbose_flag=1;
  scaling_p = vector<Probability>(T+1,Probability(0.0));//obs.length()-k+1
  //cout<<"scaling_p size: "<<scaling_p.size()<<endl;
	forward_p = vector<vector<Probability> >(nstates,vector<Probability>(T+1,Probability(0.0)));//obs.length()-k+1
	string kmer = "";

  if(verbose_flag){cout<<"//Initialize base cases t==0"<<endl;}
  if(verbose_flag){cout<<0<<'\t'<<kmer;}
  for(int s=0; s<nstates; s++)
  {
    forward_p[s][0] += start_p[s];
    scaling_p[0] += forward_p[s][0];
    if(verbose_flag){cout<<'\t'<<forward_p[s][0];}
  }
	if(verbose_flag){cout<<endl;}

  if(verbose_flag){cout<<"//Initialize base cases t==1"<<endl;}
  kmer = obs.substr(0,k);
  if(verbose_flag){cout<<1<<'\t'<<kmer;}
  // sum
	for(int s=0; s<nstates; s++)
	{
    forward_p[s][1] += forward_p[s][0]*emit_p[s][kmer];
    scaling_p[1] += forward_p[s][1];
	}
  // scale
  for(int s=0; s<nstates; s++)
  {
    //forward_p[s][1]/=scaling_p[1];
    if(verbose_flag){cout<<'\t'<<forward_p[s][1];}
  }
	if(verbose_flag){cout<<endl;}

	if(verbose_flag){cout<<"//Run forward for t>1"<<endl;}
	for(int t=2; t<=T; t++)
	{
    kmer = obs.substr(t-1,k);//t+1
		if(verbose_flag){cout<<t<<'\t'<<kmer;}
    // sum
		for(int s=0; s<nstates; s++)
		{
			for(int s0=0; s0<nstates; s0++)
			{
        forward_p[s][t] += forward_p[s0][t-1]*trans_p[s0][s]*emit_p[s][kmer];
			}
      scaling_p[t] += forward_p[s][t];
		}
    // scale
    for(int s=0; s<nstates; s++)
    {
      //forward_p[s][t]/=scaling_p[t];
      if(verbose_flag){cout<<'\t'<<forward_p[s][t];}//<<'\t'<<scaling_p[t];}
    }
		if(verbose_flag){cout<<endl;}
	}
  if(verbose_flag)
  {
    for(int t=0; t<=T; t++)//<
    {
      cout<<t<<'\t'<<scaling_p[t];
      cout<<endl;
    }
  }//verbose_flag=0;
}

void Hmm::backward(string &obs)
{
  scaling_b = vector<Probability>(T+1,Probability(0.0));//obs.length()-k+1
  //verbose_flag = 0;
  backward_p = vector<vector<Probability> >(nstates,vector<Probability>(T+1,Probability(0.0)));//init?//obs.length()-k+1
	string kmer = "";

	if(verbose_flag){cout<<"//Initialize base cases t==T"<<endl;}
	//int T = obs.length()-k;
	//kmer = obs.substr(T,k);
  kmer = obs.substr(T,k);
  //if(verbose_flag){cout<<T<<'\t'<<kmer;}
  if(verbose_flag){cout<<T<<'\t'<<kmer;}
  // end probability always unscaled 1.0
	for(int s=0; s<nstates; s++)
  {
  	//backward_p[s][T] = Probability(1.0);
  	//if(verbose_flag){cout<<'\t'<<backward_p[s][T];}
    backward_p[s][T] = Probability(1.0);
    scaling_b[T] += backward_p[s][T];
    if(verbose_flag){cout<<'\t'<<backward_p[s][T];}
  }
	if(verbose_flag){cout<<endl;}

	if(verbose_flag){cout<<"//Run backward for t<T"<<endl;}
	//for(int t=T-1; t>=0; t--)
	for(int t=T-1; t>=0; t--)//T-2
	{
    kmer = obs.substr(t,k);
		if(verbose_flag){cout<<t<<'\t'<<kmer;}
    //kmer = obs.substr(t+1,k);
		for(int s0=0; s0<nstates; s0++)
		{
			for(int s=0; s<nstates; s++)
			{
        backward_p[s0][t] += backward_p[s][t+1]*trans_p[s0][s]*emit_p[s][kmer];///scaling_p[t+1];
        scaling_b[t] += backward_p[s0][t];
			}
      if(verbose_flag){cout<<'\t'<<backward_p[s0][t];}//<<'\t'<<scaling_p[t+1];}
		}
		if(verbose_flag){cout<<endl;}
    /*//scaling could happen after summing
    for(int s0=0; s0<nstates; s0++)
    {
      backward_p[s0][t]/=scaling_p[t+1];
      if(verbose_flag){cout<<'\t'<<backward_p[s0][t];}
    }
		if(verbose_flag){cout<<endl;}*/
  }
  //verbose_flag = 0;
}

//This is the posterior probability of being in state s at time t
void Hmm::gamma(string &obs)
{
  scaling_g = vector<Probability>(T+1,Probability(0.0));//obs.length()-k+1
  gamma_p = vector<vector<Probability> >(nstates,vector<Probability>(T+1));//init?//obs.length()-k+1

	if(verbose_flag){cout<<"//Run gamma for all t"<<endl;}
	//int T = obs.length()-k+1;//forward_p[0].size();//
	for(int t=0; t<=T; t++)
	{
		if(verbose_flag){cout<<t;}
		for(int s=0; s<nstates; s++)
		{
		  //f = accumulate(forward_p[s0].begin()+1,forward_p[s0].end(),Probability(0.0),sum_prob);
  		//b = accumulate(backward_p[s0].begin()+1,backward_p[s0].end(),Probability(0.0),sum_prob);
			gamma_p[s][t] = forward_p[s][t]*backward_p[s][t];
      scaling_g[t] += gamma_p[s][t];
			if(verbose_flag){cout<<'\t'<<gamma_p[s][t];}
		}
    for(int s=0; s<nstates; s++)
    {
      gamma_p[s][t] /= scaling_g[t];
			if(verbose_flag){cout<<'\t'<<gamma_p[s][t];}
    }
		if(verbose_flag){cout<<endl;}
	}
}

void Hmm::xi(string &obs)
{
	xi_p = vector< vector< vector<Probability> > >(nstates,vector< vector<Probability> >(nstates, vector<Probability>(T)));//init?//obs.length()-k+1
  Probability scaling = Probability(0.0);
  vector<Probability> scaling_v = vector<Probability>(nstates, Probability(0.0));

	if(verbose_flag){cout<<"//Run xi for all t"<<endl;}
	string kmer = "";
	//int T = obs.length()-k+1;//forward_p[0].size();//
	for(int t=0; t<T; t++)
	{
    scaling = Probability(0.0);
    scaling_v = vector<Probability>(nstates,Probability(0.0));
		if(verbose_flag){cout<<t<<endl;}
		kmer = obs.substr(t,k);//t+1
		for(int s0=0; s0<nstates; s0++)
		{
			if(verbose_flag){cout<<s0;}
			for(int s=0; s<nstates; s++)
			{
				xi_p[s0][s][t] = forward_p[s0][t]*trans_p[s0][s]*backward_p[s][t+1]*emit_p[s][kmer]/scaling_p[T];//*scaling_p[t]*scaling_p[t+1];///scaling_p[T-1];
        scaling += xi_p[s0][s][t];
        scaling_v[s] += xi_p[s0][s][t];
				if(verbose_flag){cout<<'\t'<<xi_p[s0][s][t];}
			}
			if(verbose_flag){cout<<endl;}
		}
    //scale
		for(int s=0; s<nstates; s++)
		{
			if(verbose_flag){cout<<s;}
			for(int s0=0; s0<nstates; s0++)
			{
				//xi_p[s0][s][t] /= scaling;
  			//xi_p[s0][s][t] /= scaling_v[s];
				if(verbose_flag){cout<<'\t'<<xi_p[s0][s][t];}
			}
			if(verbose_flag){cout<<endl;}
		}
	}
}

void Hmm::update(string &obs)
{
  update_start_p();
  update_trans_p();
  update_emit_p(obs);
}

void Hmm::update_start_p()
{
	for(int s=0; s<nstates; s++)
	{
		start_p[s] = gamma_p[s][1];
	}
  start_p.normalize();
  if(verbose_flag){start_p.print();}
}

void Hmm::update_trans_p()
{
  Probability gamma = Probability(0.0);
  Probability xi = Probability(0.0);
	for(int s0=0; s0<nstates; s0++)
	{
    /*gamma = Probability(0.0);
    for(int t=0; t<gamma_p[s0].size(); t++)
    {
      gamma += gamma_p[s0][t];
    }*/
		gamma = accumulate(gamma_p[s0].begin()+1,gamma_p[s0].end()-1,Probability(0.0),sum_prob);
    // cout<<s0<<" gamma: "<<gamma<<endl;
		for(int s=0; s<nstates; s++)
		{
      /*xi = Probability(0.0);
      for(int t=0; t<gamma_p[s].size(); t++)
      {
        xi += xi_p[s0][s][t];
      }*/
  		xi = accumulate(xi_p[s0][s].begin()+1,xi_p[s0][s].end(),Probability(0.0),sum_prob);
      // cout<<"update: "<<xi<<' '<<gamma<<endl;
      trans_p[s0][s] = xi/gamma;
		}
	}
  trans_p.normalize();
  if(verbose_flag){trans_p.print();}
}

void Hmm::update_emit_p(string &obs)
{
  // Emission emit_p_new = Emission(nstates, k, emit_p.get_min());
  Emission emit_p_new = Emission(nstates, k);
  // Probability gamma_min = Probability(1.0);
  // for(int i=0; i<1; i++)//nstates
  // {
  //   gamma_min = min(gamma_min,*min_element(gamma_p[i].begin()+1,gamma_p[i].end()));
  // }
  // Emission emit_p_new = Emission(nstates, k, Probability(1.0/emit_p.get_nemissions()));
  // Emission emit_p_new = Emission(nstates, k, Probability(0.1));
  // Emission emit_p_new = Emission(nstates, k, Probability(1.0/(obs.length()-k)));
  // Emission emit_p_new = Emission(nstates, k, gamma_min);
  // Emission emit_p_new = Emission(nstates, k ,Probability(0.01));
  // cout<<1.0/(obs.length()-k)<<" "<<1.0/emit_p.get_nemissions()<<" "<<gamma_min;
	for(int s=0; s<nstates; s++)
	{
    //add gamma for each emission
		// Probability gamma = accumulate(gamma_p[s].begin()+1,gamma_p[s].end()-1,Probability(0.0),sum_prob);
    //Probability gamma_tot = Probability(0.0);
		for(int t=0; t<T; t++)//obs.length()-k+1
		{
      emit_p_new[s][obs.substr(t,k)] += gamma_p[s][t+1];///gamma;
      //gamma_tot += gamma_p[s][t];
		}
    //gamma_tot/=Probability(T);
    //cout<<" "<<gamma_tot;
    // gamma_min = *min_element(gamma_p[s].begin()+1,gamma_p[s].end());
    //emit_p_new.add(s,Probability(1.0)/gamma);
    //emit_p_new.add(s,gamma_min);
    // cout<<" "<<gamma_min;
    emit_p_new.normalize(s);

    vector<string> keys = emit_p.get_emissions();

    //add update to previous probabilities
    for(auto&& key:keys)
    {
      //emit_p[s][key] += emit_p_new[s][key];//max(emit_p_new[s][key],emit_p.get_min());
    }
    emit_p=emit_p_new;
    emit_p.normalize(s);
	}
  // cout<<endl;
  if(verbose_flag){emit_p_new.print_vertical();}
}

string Hmm::decode(string &obs)
{
  T = obs.length()-k+1;
  viterbi_p = Viterbi(nstates, k, obs);
  vector<string> path(nstates);
  vector<string> newpath(nstates);
  string kmer = "";

  if(verbose_flag){cout<<"//Initialize base cases t==0"<<endl;}
  kmer = obs.substr(0,k);
  for(int s=0; s<nstates; s++)
  {
    viterbi_p[s][0] = start_p[s]*emit_p[s][kmer];
		path[s] = to_string(s);
		if(verbose_flag){cout<<path[s].length()<<" "<<path[s]<<" "<<viterbi_p[s][0]<<" "<<endl;}
  }

	if(verbose_flag){cout<<"//Run viterbi_p for t>0"<<endl;}
	for(int t=1; t<T; t++)//obs.length()-k+1
	{
		kmer = obs.substr(t,k);
		for(int s=0; s<nstates; s++)
		{
			Probability max_p = Probability(0.0), new_p = Probability(0.0);
			int max_s = -1;
			for(int s0=0; s0<nstates; s0++)
			{
				new_p = viterbi_p[s0][t-1]*trans_p[s0][s]*emit_p[s][kmer];
        max_s = new_p>max_p?s0:max_s;
        max_p = max_p>new_p?max_p:new_p;//Probability(fmax(max_p.get_value(),new_p.get_value()),1);
			}
			viterbi_p[s][t] = max_p;
			newpath[s] = path[max_s]+to_string(s);
			if(verbose_flag){cout<<newpath[s].length()<<" "<<newpath[s]<<" "<<max_p<<endl;}
      //viterbi_p.print_vertical();
		}
		path = newpath;
	}

	if(verbose_flag){cout<<"//Search maximum"<<endl;}
  int n = fmax(0,T-1);//obs.length()-k
	Probability max_p = Probability(0.0), new_p = Probability(0.0);
	int max_s = -1;
	for(int s=0; s<nstates; s++)
  {
    new_p = viterbi_p[s][n];
    max_s = new_p>max_p?s:max_s;
    max_p = max_p>new_p?max_p:new_p;//Probability(fmax(max_p.get_value(),new_p.get_value()),1);
  }
  /*cout<<"Path: "<<path[max_s];
  cout<<"\t";
  cout<<"Prob: "<<max_p;
  cout<<endl;*/
  if(max_p.get_value()<viterbi_prob){return "";}
  set_viterbi_path(path[max_s]);
  set_viterbi_prob(max_p.get_value());
  return path[max_s];
}

string Hmm::read_fasta(string infile)
{
  string sequence = "";
  ifstream IN(infile.c_str());
  if(IN.good())
  {
    string buffer;
    while(getline(IN, buffer))
    {
      if(buffer.at(0)=='>') continue;
      string_toupper(buffer);
      sequence.append(buffer);
    }
  }
  else
  {
    printf("no such infile: %s\n",infile.c_str());
    exit(EXIT_FAILURE);
  }
  IN.close();
  return sequence;
}

void Hmm::read_emissions_fasta(int state, string infile)
{
  //map<string,double> kmer;
  //int size = 0;

	ifstream IN(infile.c_str());
	if(IN.good())
	{
		string buffer;
		while(getline(IN,buffer))
		{
			if(buffer.at(0)=='>')
			{
				continue;
			}
      string_toupper(buffer);
			for(int i=0; i<buffer.length()-k+1;i++)
			{
				//kmer[buffer.substr(i,k)]++;
  			emit_p[state][buffer.substr(i,k)]+=Probability(1.0);
        //size++;
			}
		}
	}
  else
  {
    printf("no such infile: %s\n",infile.c_str());
    exit(EXIT_FAILURE);
  }
	IN.close();
  //emit_p.print_vertical();
  emit_p.normalize(state);

  //for(auto&& emission:kmer)
  //{
  //  emit_p[state][emission.first] = Probability(emission.second/size);
  //}
  if(verbose_flag){emit_p.print_vertical();}
}

void Hmm::write_decoding(string file)
{
  auto coutbuf = cout.rdbuf();
  ofstream OUT;

  if(file!="")
  {
    file.append(".decode.txt");
    OUT.open(file.c_str());
    cout.rdbuf(OUT.rdbuf());
  }

  cout<<viterbi_prob<<endl;
  for(int i=0; i<k-1; i++)
  {
    cout<<k;
  }
  cout<<viterbi_path<<endl;

  if(file!="")
  {
    cout.rdbuf(coutbuf);
    OUT.close();
  }
}

void Hmm::write_model(string file)
{
  auto coutbuf = cout.rdbuf();
  ofstream OUT;

  if(file!="")
  {
    file.append(".model.txt");
    OUT.open(file.c_str());
    cout.rdbuf(OUT.rdbuf());
  }

  print_kmer();
  print_states();
  print_start();
  print_transitions();
  print_emissions();

  if(file!="")
  {
    cout.rdbuf(coutbuf);
    OUT.close();
  }
}

void Hmm::read_model(string file)
{
  string header = "first line";
  int state;
  double p;
  string emission;

  ifstream IN(file.c_str());
  if(IN.good())
  {
    string buffer;
    while(getline(IN,buffer))
    {
      if(buffer.at(0)=='>')
      {
        header = buffer.substr(1,buffer.length()-1);
        cout<<"header: "<<header<<endl;
        continue;
      }
      stringstream input;
      input.str(buffer);
      if(header == "kmer")
      {
        input>>k;
        cout<<"read kmer: "<<k<<endl;
        continue;
      }
      else if(header == "states")
      {
        input>>nstates;
        start_p = Start(nstates);
        trans_p = Transition(nstates);
        emit_p = Emission(nstates,k);
        cout<<"read nstates: "<<nstates<<endl;
        continue;
      }
      else if(header == "start")
      {
        input>>state;
        input>>p;
        cout<<"read start: "<<state<<" "<<p<<endl;
        start_p[state] = Probability(p);
        continue;
      }
      else if(header == "transitions")
      {
        for(int i=0; i<nstates; i++)
        {
          getline(IN,buffer);//remove state headers
          //cout<<"|"<<buffer<<"|";
          input.clear();
          input.str(buffer);
          input>>state;
          cout<<state;
          for(int j=0; j<nstates; j++)
          {
            input>>p;
            cout<<'\t'<<p;
            trans_p[i][j] = Probability(p);
          }
          cout<<endl;
        }
        continue;
      }
      else if(header == "emissions")
      {
        while(getline(IN,buffer))//remove state headers
        {
          input.clear();
          input.str(buffer);
          input>>emission;
          cout<<emission;
          for(int i=0; i<nstates; i++)
          {
            input>>p;
            cout<<'\t'<<p;
            emit_p[i][emission] = Probability(p);
          }
          cout<<endl;
        }
        header="";
        continue;
      }
      else
      {
        input>>header;
        cout<<"can not parse: "<<header<<endl;
        abort();
      }
    }
  }
}

#endif //HMM_CPP
