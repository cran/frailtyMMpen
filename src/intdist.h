#ifndef _INTDIST_H
#define _INTDIST_H

#include <Rcpp.h>

using namespace Rcpp;

struct intParams {
  
  double a;
  double b;
  double d;
  double s;
  double por;
  double mpvf;
  
  void print() {
    Rcpp::Rcout << "A: "<< a << "B: " << b << "D: "<< d << "T: "<< s << "I: " << por << "\n";
  }
  
};

Environment twd = Environment::namespace_env("mgcv");
Function TW = twd["ldTweedie"];

double logN1int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = 1/(x*std::sqrt(2*M_PI*k.s))*std::exp(-std::pow(std::log(x), 2.0)/(2*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x);
  return ret;
  
}

double logN2int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = 1/(x*std::sqrt(2*M_PI*k.s))*std::exp(-std::pow(std::log(x), 2.0)/(2*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x)*x/k.por;
  return ret;
  
}

double logN3int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::pow(std::log(x), 2.0)/(x*std::sqrt(2*M_PI*k.s))*std::exp(-std::pow(std::log(x), 2.0)/(2*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x)/k.por;
  return ret;
  
}

double InvG1int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = 1/(std::sqrt(2*M_PI*k.s)*std::pow(x, 1.5))*std::exp(-std::pow(x-1, 2.0)/(2*x*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x);
  return ret;
  
}

double InvG2int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = 1/(std::sqrt(2*M_PI*k.s)*std::pow(x, 1.5))*std::exp(-std::pow(x-1, 2.0)/(2*x*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x)*x/k.por;
  return ret;
  
}

double InvG3int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::pow(x-1, 2.0)/x/(std::sqrt(2*M_PI*k.s)*std::pow(x, 1.5))*std::exp(-std::pow(x-1, 2.0)/(2*x*k.s))*std::pow(x, k.d)*k.b*std::exp(-k.a*x)/k.por;
  return ret;
  
}


double PVF1int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::exp(*REAL(TW(x, 1, k.mpvf, k.s)));
  ret *= std::pow(x, k.d)*k.b*std::exp(-k.a*x);
  return ret;
  
}

double PVF2int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::exp(*REAL(TW(x, 1, k.mpvf, k.s)));
  ret *= std::pow(x, k.d)*k.b*std::exp(-k.a*x)*x/k.por;
  return ret;
  
}

double PVF3int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::exp(*REAL(TW(x, 1, k.mpvf, k.s)));
  ret *= *(REAL(TW(x, 1, k.mpvf, k.s))+1);
  ret *= std::pow(x, k.d)*k.b*std::exp(-k.a*x)/k.por;
  return ret;
  
}

double PVF4int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::exp(*REAL(TW(x, 1, k.mpvf, k.s)));
  ret *= *(REAL(TW(x, 1, k.mpvf, k.s))+2);
  ret *= std::pow(x, k.d)*k.b*std::exp(-k.a*x)/k.por;
  return ret;
  
}

double PVF5int(double x, void * params) {
  
  intParams k = *(intParams *) params;
  double ret = std::exp(*REAL(TW(x, 1, k.mpvf, k.s)));
  ret *= *REAL(TW(x, 1, k.mpvf, k.s));
  ret *= std::pow(x, k.d)*k.b*std::exp(-k.a*x)/k.por;
  return ret;
  
}

#endif
