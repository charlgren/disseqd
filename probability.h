#ifndef PROBABILITY_H
#define PROBABILITY_H

#include <cmath>
#include <iostream>
#include "functions.cpp"

using namespace std;

class Probability
{
public:
//-------------------------------- constructors --------------------------------
  Probability(){value=probability(1.0);}
  Probability(double v, double logged=0){value=logged?v:probability(v);}//check(v)
//--------------------------------- operators ----------------------------------
  friend bool operator==(const Probability& lhs, const Probability& rhs){return lhs.value == rhs.value;}
  friend bool operator<(const Probability& lhs, const Probability& rhs){return lhs.value < rhs.value;}
  friend bool operator>(const Probability& lhs, const Probability& rhs){return lhs.value > rhs.value;}
  Probability& operator=(const Probability& rhs){value = rhs.value; return *this;}
  friend ostream& operator<<(ostream& out, const Probability& rhs);
  #ifdef LOGSCALE //log probabilities
    Probability operator+(const Probability& rhs){return Probability(sum_log_prob(get_value(),rhs.get_value()),1);}
    Probability operator-(const Probability& rhs){return Probability(sum_log_prob(get_value(),-rhs.get_value()),1);}
    Probability operator*(const Probability& rhs){return Probability(get_value()+rhs.get_value(),1);}
    Probability operator/(const Probability& rhs){return Probability(get_value()-rhs.get_value(),1);}
    Probability operator+=(const Probability& rhs){return set_value(sum_log_prob(get_value(),rhs.get_value()));return *this;}
    Probability operator-=(const Probability& rhs){return set_value(sum_log_prob(get_value(),-rhs.get_value()));return *this;}
    Probability operator*=(const Probability& rhs){return set_value(get_value()+rhs.get_value());return *this;}
    Probability operator/=(const Probability& rhs){return set_value(get_value()-rhs.get_value());return *this;}
  #else //std probabilities
    Probability operator+(const Probability& rhs){return Probability(sum_std_prob(get_value(),rhs.get_value()));}
    Probability operator-(const Probability& rhs){return Probability(sum_std_prob(get_value(),-rhs.get_value()));}
    Probability operator*(const Probability& rhs){return Probability(get_value()*rhs.get_value());}
    Probability operator/(const Probability& rhs){return Probability(get_value()/rhs.get_value());}
    Probability operator+=(const Probability& rhs){return set_value(sum_std_prob(get_value(),rhs.get_value()));return *this;}
    Probability operator-=(const Probability& rhs){return set_value(sum_std_prob(get_value(),-rhs.get_value()));return *this;}
    Probability operator*=(const Probability& rhs){return set_value(get_value()*rhs.get_value());return *this;}
    Probability operator/=(const Probability& rhs){return set_value(get_value()/rhs.get_value());return *this;}
  #endif //LOGSCALE
//--------------------------------- functions ----------------------------------
  double get_value() const {return value;}
  double set_value(double v){value=v; return value;}//check(v)
  //void print(){cout<<value;}

private:
//--------------------------------- variables ----------------------------------
  double value;
};

#ifdef LOGSCALE
  //log_probabilities
  static Probability sum_prob(Probability a, Probability b){return Probability(sum_log_prob(a.get_value(),b.get_value()),1);}
  ostream& operator<<(ostream& out, const Probability& rhs){out<<exp(rhs.get_value());return out;}
#else
  //std_probabilities
  static Probability sum_prob(Probability a, Probability b){return Probability(sum_std_prob(a.get_value(),b.get_value()));}
  ostream& operator<<(ostream& out, const Probability& rhs){out<<rhs.get_value();return out;}
#endif //LOGSCALE

#endif //PROBABILITY_H
