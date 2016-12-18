#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cctype>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

//DBL_EPS: 2.22045e-16 -36.0437
//DBL_MIN: 2.22507e-308 -708.396
//DBL_MAX: 1.79769e+308 709.783

//LDBL_EPS: 1.0842e-19 -43.6683
//LDBL_MIN: 3.3621e-4932 -11355.1
//LDBL_MAX: 1.18973e+4932 11356.5

static double log1pexp(double x){return x<log(DBL_MIN)? 0.0: log1p(exp(x));}

static double sum_log_prob(double a, double b){return a>b? a+log1pexp(b-a):  b+log1pexp(a-b);}
static double sum_std_prob(double a, double b){return a+b;}

static double std_prob(double x){return x<DBL_MIN? 0.0: x;}
static double log_prob(double x){return x<DBL_MIN? -DBL_MAX: log(x);}

#ifdef LOGSCALE
  //log probabilities
  static double probability(double x){return log_prob(x);}
#else
  //std probabilities
  static double probability(double x){return std_prob(x);}
#endif

vector<string> split(string str, char delimiter)
{
  vector<string> internal;
  stringstream ss(str);
  string tok;

  while(getline(ss, tok, delimiter))
  {
    internal.push_back(tok);
  }

  return internal;
}

void string_toupper(string &s)
{
  std::transform(s.begin(),s.end(),s.begin(),(int(*)(int))std::toupper);
}

#endif
