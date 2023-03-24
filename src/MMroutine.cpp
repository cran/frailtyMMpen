#include <RcppGSL.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <string>
#include <vector>
#include <algorithm>
#include "intdist.h"

using namespace Rcpp;

struct Reorder {

  double a;
  double b;
  int id;

};

void computeLAM(NumericVector& LAM, const NumericVector& lambda, const NumericVector& y, const int& N, int mode) {

  std::vector<Reorder> LaR(N);
  std::vector<double> cumLam(N);

  for (int i = 0; i < N; i++) {
    LaR[i].a = lambda[i];
    LaR[i].b = y[i];
    LaR[i].id = i;
  }

  if (mode == 0) {
    std::sort(LaR.begin(), LaR.end(), [](const Reorder& z1, const Reorder& z2){return(z1.b < z2.b);});
  }

  if (mode == 1) {
    std::sort(LaR.begin(), LaR.end(), [](const Reorder& z1, const Reorder& z2){return(z1.b > z2.b);});
  }

  for (int i = 0; i < N; i++) {
    cumLam[i] = LaR[i].a;
  }

  std::partial_sum(cumLam.begin(), cumLam.end(), cumLam.begin(), std::plus<double>());

  for (int i = 0; i < N; i++) {
    LAM[LaR[i].id] = cumLam[i];
  }

}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


// [[Rcpp::export]]
List MMCL(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0,
          const double& tht0, int frailty, int penalty, double gam, double tune, const NumericVector& id, int N, int a, int p, double power) {
  
  NumericVector coef = clone(coef0);
  NumericVector lambda = clone(lambda0);
  double tht = tht0;
  
  NumericVector La(N);
  
  NumericVector YpreExp(N, 0.0);
  intParams kint;
  
  NumericVector AA(N);
  NumericVector BB(N);
  
  NumericVector A(a, 0.0);
  NumericVector B(a, 1.0);
  NumericVector D(a, 0.0);
  
  NumericVector AT(a, 0.0);
  NumericVector DT(a, 0.0);
  
  NumericVector int1(a, 0.0);
  NumericVector int2(a, 0.0);
  NumericVector int3(a, 0.0);
  
  NumericVector ME(N, 0.0);
  NumericVector E0(N, 0.0);
  NumericVector SUM0(N, 0.0);
  
  NumericVector AVEX(N, 0.0);
  NumericVector E1(N, 0.0);
  NumericVector D11(N, 0.0);
  NumericVector SUM1(N, 0.0);
  NumericVector E2(N, 0.0);
  NumericVector SUM2(N, 0.0);
  NumericVector D22(N, 0.0);
  
  int tempID(0);
  double D1(0.0), D2(0.0);
  
  computeLAM(La, lambda, y, N, 0);
  
  for (int i = 0; i < p; i++) {
    YpreExp += X[Range(i*N, (i+1)*N-1)] * coef[i];
  }
  
  YpreExp = exp(YpreExp);
  
  AA = La * YpreExp;
  BB = lambda * YpreExp;
  
  for (int j = 0; j < N; j++) {
    tempID = id[j];
    A[tempID] += AA[j];
    if (d[j] == 1) {
      B[tempID] *= BB[j];
      D[tempID]++;
    }
  }
  
  if (frailty == 1 || frailty == 2) {
    
    gsl_set_error_handler_off();
    int status;
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1, F2, F3;
    double result(0.0), error;
    
    if (frailty == 1) {
      F1.function = &logN1int;
      F2.function = &logN2int;
      F3.function = &logN3int;
    }
    
    if (frailty == 2) {
      F1.function = &InvG1int;
      F2.function = &InvG2int;
      F3.function = &InvG3int;
    }
    
    F1.params = &kint;
    F2.params = &kint;
    F3.params = &kint;
    kint.s = tht;
    
    for (int i = 0; i < a; i++) {
      
      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];
      
      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        
        status = gsl_integration_qags (&F1, 0.001, 10, 0, 1e-7, 1000, w, &result, &error);
        Rcout << "F1: Approximated by Interval\n";
        
      }
      
      kint.por = result;
      int1[i] = result;
      
      status = gsl_integration_qagiu (&F2, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        Rcout << "F2: 0\n";
        return(List::create(_["error"] = 1));
      }
      
      int2[i] = result;
      
      status = gsl_integration_qagiu (&F3, 0, 0, 1e-7, 1000, w, &result, &error);
      int3[i] = result;
      
    }
    
    tht = std::accumulate(int3.begin(), int3.end(), 0.0) / a;
    
    gsl_integration_workspace_free (w);
    
  }
  
  if (frailty == 3) {
    
    NumericVector int4(a, 0.0);
    
    gsl_set_error_handler_off();
    
    int status;
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1, F2, F3, F4;
    double result(0.0), error;
    
    F1.function = &PVF1int;
    F2.function = &PVF2int;
    F3.function = &PVF3int;
    F4.function = &PVF4int;
    
    F1.params = &kint;
    F2.params = &kint;
    F3.params = &kint;
    F4.params = &kint;
    kint.s = tht;
    kint.mpvf = power;
    
    for (int i = 0; i < a; i++) {
      
      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];
      
      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);
      
      kint.por = result;
      int1[i] = result;
      
      status = gsl_integration_qagiu (&F2, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        Rcout << "F2: 0\n";
        return(List::create(_["error"] = 1));
      }
      
      int2[i] = result;
      
      status = gsl_integration_qagiu (&F3, 0, 0, 1e-7, 1000, w, &result, &error);
      int3[i] = result;
      
      status = gsl_integration_qagiu (&F4, 0, 0, 1e-7, 1000, w, &result, &error);
      int4[i] = result;
      
    }
    
    double thtPVF = tht - std::accumulate(int3.begin(), int3.end(), 0.0) / std::accumulate(int4.begin(), int4.end(), 0.0);
    
    gsl_integration_workspace_free (w);
    
    if (thtPVF > 0) {
      tht = thtPVF;
    }
    
  }
  
  if (frailty == 0) {
    
    AT = 1/tht + A;
    DT = 1/tht + D;
    int2 = DT/AT;
    
    double tht1(0.0), Q1(0.0), Q21(0.0), Q22(0.0);
    
    Q1 += a*(gsl_sf_psi(1/tht) + std::log(tht) - 1);
    for (int i = 0; i < a; i++) {
      Q1 += int2[i] + std::log(AT[i]) - gsl_sf_psi(DT[i]);
    }
    Q1 /= std::pow(tht, 2);
    
    Q21 += a*(3 - 2*gsl_sf_psi(1/tht) - 2*std::log(tht));
    for (int i = 0; i < a; i++) {
      Q21 += 2*(gsl_sf_psi(DT[i]) - int2[i] - std::log(AT[i]));
    }
    Q21 /= std::pow(tht, 3);
    
    Q22 = - a*gsl_sf_psi_1(1/tht) / std::pow(tht, 4);
    
    tht1 = tht - Q1 / (Q21 + Q22);
    if(tht1 > 0) {
      tht = tht1 ;
    }
    
  }
  
  // Update lambda Variables
  
  for (int j = 0; j < N; j++) {
    tempID = id[j];
    ME[j] = int2[tempID];
  }
  
  E0 = ME * YpreExp;
  computeLAM(SUM0, E0, y, N, 1);
  lambda = d / SUM0;
  
  // Update coefficients
  
  for (int i = 0; i < p; i++) {
    AVEX += abs(X[Range(i*N, (i+1)*N-1)]);
  } 
  
  for (int i = 0; i < p; i++) {
    
    E1 = ME * X[Range(i*N, (i+1)*N-1)] * YpreExp;
    E2 = ME * AVEX *X[Range(i*N, (i+1)*N-1)] * YpreExp;
    
    computeLAM(SUM1, E1, y, N, 1);
    computeLAM(SUM2, E2, y, N, 1);
    
    D11 = d * X[Range(i*N, (i+1)*N-1)] - d * SUM1 / SUM0;
    D22 = - d * SUM2 / SUM0;
    
    D1 = std::accumulate(D11.begin(), D11.end(), 0.0);
    D2 = std::accumulate(D22.begin(), D22.end(), 0.0);
    
    if (penalty != 0) {
      
      if (penalty == 1) {
        D1 -= N*sgn(coef[i])*tune;
        D2 -= N*tune/std::abs(coef[i]);
      }
      
      if (penalty == 2) {
        D1 -= N*sgn(coef[i])*(tune - std::abs(coef[i])/gam)*(std::abs(coef[i]) <= gam*tune);
        D2 -= N*(tune - std::abs(coef[i])/gam)*(std::abs(coef[i]) <= gam*tune)/std::abs(coef[i]);
      }
      
      if (penalty == 3) {
        D1 -= N*sgn(coef[i])*(tune*(std::abs(coef[i]) <= tune) + std::max(0.0, gam*tune - std::abs(coef[i]))*(std::abs(coef[i]) > tune)/(gam - 1));
        D2 -= N*(tune*(std::abs(coef[i]) <= tune) + std::max(0.0, gam*tune - std::abs(coef[i]))*(std::abs(coef[i]) > tune)/(gam - 1))/std::abs(coef[i]);
      }
    }
    
    coef[i] -= D1/D2;
    
  }
  
  List ret = List::create(_["coef"] = coef, _["est.tht"] = tht, _["lambda"] = lambda, _["error"] = 0, _["Ar"] = int2);
  return ret;
}

// [[Rcpp::export]]
List MMME(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0,
          const double& tht0, int frailty, int penalty, double gam, double tune, int N, int n, int p, double power) {
  
  NumericVector coef = clone(coef0);
  NumericVector lambda = clone(lambda0);
  NumericVector dn = clone(d);
  double tht = tht0;
  
  int ne = N / n;
  
  NumericVector YpreExp(N);
  NumericVector La(N);
  
  intParams kint;
  
  NumericVector A(n, 0.0);
  NumericVector B(n, 1.0);
  NumericVector D(n, 0.0);
  
  NumericVector AT(n, 0.0);
  NumericVector DT(n, 0.0);
  
  NumericVector int1(n, 0.0);
  NumericVector int2(n, 0.0);
  NumericVector int3(n, 0.0);
  
  NumericVector E0n(N);
  NumericVector SUM0n(N);
  NumericVector AVEXn(N);
  
  for (int i = 0; i < ne; i++) {
    
    Range RN = Range(i*n, (i+1)*n-1);
    NumericVector LaTemp(n);
    computeLAM(LaTemp, lambda[RN], y[RN], n, 0);
    std::copy(LaTemp.begin(), LaTemp.end(), La.begin()+i*n);

    for (int j = 0; j < p; j++) {
      YpreExp[RN] += X[Range(i*n+j*N, (i+1)*n-1+j*N)] * coef[j];
    }
    
  }
  YpreExp = exp(YpreExp);
  
  for (int i = 0; i < ne; i++) {
    
    Range RN = Range(i*n, (i+1)*n-1);
    A += La[RN] * YpreExp[RN];
    for (int j = 0; j < n; j++) {
      if (dn[i*n + j] == 1) {
        B[j] *= lambda[i*n + j] * YpreExp[i*n + j];
      }
    }
    D += dn[RN];
  }
  
  if (frailty == 1 || frailty == 2) {

    gsl_set_error_handler_off();
    int status;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1, F2, F3;
    double result(0.0), error;

    if (frailty == 1) {
      F1.function = &logN1int;
      F2.function = &logN2int;
      F3.function = &logN3int;
    }

    if (frailty == 2) {
      F1.function = &InvG1int;
      F2.function = &InvG2int;
      F3.function = &InvG3int;
    }

    F1.params = &kint;
    F2.params = &kint;
    F3.params = &kint;
    kint.s = tht;

    for (int i = 0; i < n; i++) {

      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];

      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);

      if (status) {

        status = gsl_integration_qags (&F1, 0.001, 10, 0, 1e-7, 1000, w, &result, &error);
        Rcout << "F1: Approximated by Interval\n";

      }

      kint.por = result;
      int1[i] = result;

      status = gsl_integration_qagiu (&F2, 0, 0, 1e-7, 1000, w, &result, &error);

      if (status) {
        Rcout << "F2: 0\n";
        return(List::create(_["error"] = 1));
      }

      int2[i] = result;

      status = gsl_integration_qagiu (&F3, 0, 0, 1e-7, 1000, w, &result, &error);
      int3[i] = result;

    }

    tht = std::accumulate(int3.begin(), int3.end(), 0.0) / n;

    gsl_integration_workspace_free (w);

  }

  if (frailty == 3) {

    NumericVector int4(n, 0.0);

    gsl_set_error_handler_off();

    int status;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1, F2, F3, F4;
    double result(0.0), error;

    F1.function = &PVF1int;
    F2.function = &PVF2int;
    F3.function = &PVF3int;
    F4.function = &PVF4int;

    F1.params = &kint;
    F2.params = &kint;
    F3.params = &kint;
    F4.params = &kint;
    kint.s = tht;
    kint.mpvf = power;

    for (int i = 0; i < n; i++) {

      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];

      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);

      kint.por = result;
      int1[i] = result;

      status = gsl_integration_qagiu (&F2, 0, 0, 1e-7, 1000, w, &result, &error);

      if (status) {
        Rcout << "F2: 0\n";
        return(List::create(_["error"] = 1));
      }

      int2[i] = result;

      status = gsl_integration_qagiu (&F3, 0, 0, 1e-7, 1000, w, &result, &error);
      int3[i] = result;

      status = gsl_integration_qagiu (&F4, 0, 0, 1e-7, 1000, w, &result, &error);
      int4[i] = result;

    }

    double thtPVF = tht - std::accumulate(int3.begin(), int3.end(), 0.0) / std::accumulate(int4.begin(), int4.end(), 0.0);

    gsl_integration_workspace_free (w);

    if (thtPVF > 0) {
      tht = thtPVF;
    }

  }

  if (frailty == 0) {

    AT = 1/tht + A;
    DT = 1/tht + D;
    int2 = DT/AT;

    double tht1(0.0), Q1(0.0), Q21(0.0), Q22(0.0);

    Q1 += n*(gsl_sf_psi(1/tht) + std::log(tht) - 1);
    for (int i = 0; i < n; i++) {
      Q1 += int2[i] + std::log(AT[i]) - gsl_sf_psi(DT[i]);
    }
    Q1 /= std::pow(tht, 2);

    Q21 += n*(3 - 2*gsl_sf_psi(1/tht) - 2*std::log(tht));
    for (int i = 0; i < n; i++) {
      Q21 += 2*(gsl_sf_psi(DT[i]) - int2[i] - std::log(AT[i]));
    }
    Q21 /= std::pow(tht, 3);

    Q22 = - n*gsl_sf_psi_1(1/tht) / std::pow(tht, 4);

    tht1 = tht - Q1 / (Q21 + Q22);
    if(tht1 > 0) {
      tht = tht1 ;
    }

  }

  // Update lambda Variables

  for (int i = 0; i < ne; i++) {
    
    Range RN = Range(i*n, (i+1)*n-1);
    NumericVector LaTemp(n);
    
    E0n[RN] = int2 * YpreExp[RN];
    computeLAM(LaTemp, E0n[RN], y[RN], n, 1);
    std::copy(LaTemp.begin(), LaTemp.end(), SUM0n.begin()+i*n);
    lambda[RN] = dn[RN] / SUM0n[RN];
    
    for (int j = 0; j < p; j++) {
      AVEXn[RN] += abs(X[Range(i*n+j*N, (i+1)*n-1+j*N)]);
    }
    
  }

  // Update coefficients

  for (int i = 0; i < p; i++) {

    NumericVector E1(n, 0.0);
    NumericVector E2(n, 0.0);
    
    double D1(0.0), D2(0.0);
    
    for (int j = 0; j < ne; j++) {
      
      Range RN = Range(j*n, (j+1)*n-1);
      E1 += La[RN] * X[Range(j*n+i*N, (j+1)*n-1+i*N)] * YpreExp[RN];
      D1 += sum(dn[RN]*X[Range(j*n+i*N, (j+1)*n-1+i*N)]);
      E2 += La[RN] * abs(X[Range(j*n+i*N, (j+1)*n-1+i*N)]) * AVEXn[RN] * YpreExp[RN];
    }

    D1 -= sum(int2*E1);

    if (frailty != 0) {
      D2 = -2*sum(int2*E2);
    }
    else {
      D2 = -2*sum((D + 2/tht)*E2/AT);
    }

    if (penalty != 0) {

      if (penalty == 1) {
        D1 -= n*sgn(coef[i])*tune;
        D2 -= n*tune/std::abs(coef[i]);
      }

      if (penalty == 2) {
        D1 -= n*sgn(coef[i])*(tune - std::abs(coef[i])/gam)*(std::abs(coef[i]) <= gam*tune);
        D2 -= n*(tune - std::abs(coef[i])/gam)*(std::abs(coef[i]) <= gam*tune)/std::abs(coef[i]);
      }

      if (penalty == 3) {
        D1 -= n*sgn(coef[i])*(tune*(std::abs(coef[i]) <= tune) + std::max(0.0, gam*tune - std::abs(coef[i]))*(std::abs(coef[i]) > tune)/(gam - 1));
        D2 -= n*(tune*(std::abs(coef[i]) <= tune) + std::max(0.0, gam*tune - std::abs(coef[i]))*(std::abs(coef[i]) > tune)/(gam - 1))/std::abs(coef[i]);
      }
    }

    coef[i] -= D1/D2;

  }
  
  List ret = List::create(_["coef"] = coef, _["est.tht"] = tht, _["lambda"] = lambda, _["error"] = 0, _["Ar"] = int2);
  return ret;
}



// [[Rcpp::export]]
List MMRE(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0,
          const double& tht0, int frailty, int penalty, double gam, double tune, const NumericVector& id, int N, int a, int p, double power) {
  
  NumericVector coef = clone(coef0);
  NumericVector lambda = clone(lambda0);
  double tht = tht0;
  
  NumericVector La(N);
  
  NumericVector XM(a*p*N, 0.0);
  NumericVector XMexp(a*N, 0.0);
  NumericVector YpreExp(N, 0.0);
  NumericVector Yend(a, 0.0);
  intParams kint;
  
  NumericVector A(a, 0.0);
  NumericVector B(a, 1.0);
  NumericVector D(a, 0.0);
  
  NumericVector AT(a, 0.0);
  NumericVector DT(a, 0.0);
  
  NumericVector int1(a, 0.0);
  NumericVector int2(a, 0.0);
  NumericVector int3(a, 0.0);
  
  NumericVector SUM0(N, 0.0);
  
  NumericVector AVEX(a*N, 0.0);
  NumericVector SUM1(N, 0.0);
  NumericVector SUM2(N, 0.0);
  
  int tempID(0);
  double D1(0.0), D2(0.0);
  
  for (int i = 0; i < N; i++) {
    tempID = id[i];
    if (i == N-1 || tempID < id[i+1]) {
      Yend[tempID] = y[i];
    }
  }
  
  computeLAM(La, lambda, y, N, 0);
  
  for (int i = 0; i < p; i++) {
    YpreExp += X[Range(i*N, (i+1)*N-1)] * coef[i];
  }
  
  YpreExp = exp(YpreExp);
  
  for (int i = N - 1; i >= 0; i--) {
    tempID = id[i];
    double ct = y[i];
    if (i == N - 1 || tempID < id[i+1]) {
      for (int j = 0; j < N; j++) {
        for (int z = 0; z < p; z++) {
          XM[j*a*p + z*a + tempID] = X[z*N + i];
        }
      }
    } else {
      for (int j = 0; j < N; j++) {
        if (ct >= y[j]) {
          for (int z = 0; z < p; z++) {
            XM[j*a*p + z*a + tempID] = X[z*N + i];
          }
        }
      }
    }
  }
  
  for (int i = 0; i < N; i++) {
    for (int z = 0; z < p; z++) {
      XMexp[Range(i*a, (i+1)*a-1)] += XM[Range(i*a*p+z*a, i*a*p+(z+1)*a-1)] * coef[z];
    }
  }
  XMexp = exp(XMexp);
  
  double L0 = 0;
  for (int i = 0; i < N; i++) {
    tempID = id[i];
    if (i == 0 || tempID > id[i-1]) {
      L0 = 0;
    }
    A[tempID] += (La[i] - L0) * YpreExp[i];
    L0 = La[i];
  }
  
  for (int j = 0; j < N; j++) {
    tempID = id[j];
    if (d[j] == 1) {
      B[tempID] *= lambda[j] * YpreExp[j];
      D[tempID]++;
    }
  }
  
  if (frailty == 1 || frailty == 2) {
    
    gsl_set_error_handler_off();
    int status;
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1, F2, F3;
    double result(0.0), error;
    
    if (frailty == 1) {
      F1.function = &logN1int;
      F2.function = &logN2int;
      F3.function = &logN3int;
    }
    
    if (frailty == 2) {
      F1.function = &InvG1int;
      F2.function = &InvG2int;
      F3.function = &InvG3int;
    }
    
    F1.params = &kint;
    F2.params = &kint;
    F3.params = &kint;
    kint.s = tht;
    
    for (int i = 0; i < a; i++) {
      
      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];
      
      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        
        status = gsl_integration_qags (&F1, 0.001, 10, 0, 1e-7, 1000, w, &result, &error);
        Rcout << "F1: Approximated by Interval\n";
        
      }
      
      kint.por = result;
      int1[i] = result;
      
      status = gsl_integration_qagiu (&F2, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        Rcout << "F2: 0\n";
        return(List::create(_["error"] = 1));
      }
      
      int2[i] = result;
      
      status = gsl_integration_qagiu (&F3, 0, 0, 1e-7, 1000, w, &result, &error);
      int3[i] = result;
      
    }
    
    tht = std::accumulate(int3.begin(), int3.end(), 0.0) / a;
    
    gsl_integration_workspace_free (w);
    
  }
  
  if (frailty == 3) {
    
    NumericVector int4(a, 0.0);
    
    gsl_set_error_handler_off();
    
    int status;
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1, F2, F3, F4;
    double result(0.0), error;
    
    F1.function = &PVF1int;
    F2.function = &PVF2int;
    F3.function = &PVF3int;
    F4.function = &PVF4int;
    
    F1.params = &kint;
    F2.params = &kint;
    F3.params = &kint;
    F4.params = &kint;
    kint.s = tht;
    kint.mpvf = power;
    
    for (int i = 0; i < a; i++) {
      
      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];
      
      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);
      
      kint.por = result;
      int1[i] = result;
      
      status = gsl_integration_qagiu (&F2, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        Rcout << "F2: 0\n";
        return(List::create(_["error"] = 1));
      }
      
      int2[i] = result;
      
      status = gsl_integration_qagiu (&F3, 0, 0, 1e-7, 1000, w, &result, &error);
      int3[i] = result;
      
      status = gsl_integration_qagiu (&F4, 0, 0, 1e-7, 1000, w, &result, &error);
      int4[i] = result;
      
    }
    
    double thtPVF = tht - std::accumulate(int3.begin(), int3.end(), 0.0) / std::accumulate(int4.begin(), int4.end(), 0.0);
    
    gsl_integration_workspace_free (w);
    
    if (thtPVF > 0) {
      tht = thtPVF;
    }
    
  }
  
  if (frailty == 0) {
    
    AT = 1/tht + A;
    DT = 1/tht + D;
    int2 = DT/AT;
    
    double tht1(0.0), Q1(0.0), Q21(0.0), Q22(0.0);
    
    Q1 += a*(gsl_sf_psi(1/tht) + std::log(tht) - 1);
    for (int i = 0; i < a; i++) {
      Q1 += int2[i] + std::log(AT[i]) - gsl_sf_psi(DT[i]);
    }
    Q1 /= std::pow(tht, 2);
    
    Q21 += a*(3 - 2*gsl_sf_psi(1/tht) - 2*std::log(tht));
    for (int i = 0; i < a; i++) {
      Q21 += 2*(gsl_sf_psi(DT[i]) - int2[i] - std::log(AT[i]));
    }
    Q21 /= std::pow(tht, 3);
    
    Q22 = - a*gsl_sf_psi_1(1/tht) / std::pow(tht, 4);
    
    tht1 = tht - Q1 / (Q21 + Q22);
    if(tht1 > 0) {
      tht = tht1 ;
    }
    
  }
  
  // Update lambda Variables
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < a; j++) {
      SUM0[i] += (Yend[j] >= y[i]) * int2[j] * XMexp[i*a+j];
    }
  }
  
  lambda = d / SUM0;
  
  // Update coefficients
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < a; j++) {
      for (int z = 0; z < p; z++) {
        AVEX[i*a + j] += std::abs(XM[i*a*p + z*a + j]);
      } 
    }
  }
  
  for (int i = 0; i < p; i++) {
    
    for (int z = 0; z < N; z++) {
      SUM1[z] = 0;
      SUM2[z] = 0;
    }
    
    for (int z = 0; z < N; z++) {
      for (int j = 0; j < a; j++) {
        SUM1[z] += (Yend[j] >= y[z]) * int2[j] * XMexp[z*a+j] * XM[z*a*p+i*a+j];
        SUM2[z] += (Yend[j] >= y[z]) * int2[j] * XMexp[z*a+j] * std::abs(XM[z*a*p+i*a+j]) * AVEX[z*a+j];
      }
    }
    
    D1 = 0;
    D2 = 0;
    for (int z = 0; z < N; z++) {
      D1 += d[z] * X[i*N + z] - d[z] * SUM1[z] / SUM0[z];
      D2 -= d[z] * SUM2[z] / SUM0[z];
    }
    
    if (penalty != 0) {
      
      if (penalty == 1) {
        D1 -= N*sgn(coef[i])*tune;
        D2 -= N*tune/std::abs(coef[i]);
      }
      
      if (penalty == 2) {
        D1 -= N*sgn(coef[i])*(tune - std::abs(coef[i])/gam)*(std::abs(coef[i]) <= gam*tune);
        D2 -= N*(tune - std::abs(coef[i])/gam)*(std::abs(coef[i]) <= gam*tune)/std::abs(coef[i]);
      }
      
      if (penalty == 3) {
        D1 -= N*sgn(coef[i])*(tune*(std::abs(coef[i]) <= tune) + std::max(0.0, gam*tune - std::abs(coef[i]))*(std::abs(coef[i]) > tune)/(gam - 1));
        D2 -= N*(tune*(std::abs(coef[i]) <= tune) + std::max(0.0, gam*tune - std::abs(coef[i]))*(std::abs(coef[i]) > tune)/(gam - 1))/std::abs(coef[i]);
      }
    }
    
    coef[i] -= D1/D2;
    
  }
  
  List ret = List::create(_["coef"] = coef, _["est.tht"] = tht, _["lambda"] = lambda, _["error"] = 0, _["Ar"] = int2);
  return ret;
}


// [[Rcpp::export]]
double LogLikCL(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0,
                const double& tht0, int frailty, const NumericVector& id, int N, int a, int p, double power) {
  

  NumericVector coef = clone(coef0);
  NumericVector lambda = clone(lambda0);
  double tht = tht0;
  
  NumericVector La(N);
  
  NumericVector Ypre(N, 0.0);
  NumericVector YpreExp(N, 0.0);
  intParams kint;
  
  NumericVector AA(N);
  NumericVector BB(N);
  
  NumericVector A(a, 0.0);
  NumericVector B(a, 1.0);
  NumericVector D(a, 0.0);
  
  NumericVector AT(a, 0.0);
  NumericVector DT(a, 0.0);
  
  NumericVector int1(a, 0.0);
  
  int tempID(0);
  
  computeLAM(La, lambda, y, N, 0);
  
  for (int i = 0; i < p; i++) {
    Ypre += X[Range(i*N, (i+1)*N-1)] * coef[i];
  }
  
  YpreExp = exp(Ypre);
  
  AA = La * YpreExp;
  BB = lambda * YpreExp;
  
  for (int j = 0; j < N; j++) {
    tempID = id[j];
    A[tempID] += AA[j];
    if (d[j] == 1) {
      B[tempID] *= BB[j];
      D[tempID]++;
    }
  }
  
  if (frailty == 1 || frailty == 2) {
    
    gsl_set_error_handler_off();
    int status;
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1;
    double result(0.0), error;
    
    if (frailty == 1) {
      F1.function = &logN1int;
    }
    
    if (frailty == 2) {
      F1.function = &InvG1int;
    }
    
    if (frailty == 3) {
      F1.function = &PVF1int;
      kint.mpvf = power;
    }
    
    F1.params = &kint;
    kint.s = tht;
    
    for (int i = 0; i < a; i++) {
      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];
      
      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        
        status = gsl_integration_qags (&F1, 0.001, 10, 0, 1e-7, 1000, w, &result, &error);
        Rcout << "F1: Approximated by Interval\n";
        
      }
      
      int1[i] = result;
      
    }
    
    gsl_integration_workspace_free (w);
    
    return(sum(log(int1)));
  }
  
  if (frailty == 0) {
    AT = 1/tht + A;
    DT = 1/tht + D;
    
    double l(0.0);
    l -= a*(lgamma(1/tht) + std::log(tht)/tht);
    l += sum(lgamma(DT)) - sum(DT*log(AT));
    l += sum(d * Ypre);
    
    for (int i = 0; i < N; i++) {
      if (lambda[i] > 0) {
        l += std::log(lambda[i]);
      }
    }
    
    return l;
    
  } 
  
  return 0.0;
}


// [[Rcpp::export]]
double LogLikME(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0,
                const double& tht0, int frailty, int N, int n, int p, double power) {
  
  
  NumericVector coef = clone(coef0);
  NumericVector lambda = clone(lambda0);
  NumericVector dn = clone(d);
  double tht = tht0;
  
  int ne = N / n;
  
  NumericVector Ypre(N);
  NumericVector YpreExp(N);
  NumericVector La(N);
  
  intParams kint;
  
  NumericVector A(n, 0.0);
  NumericVector B(n, 1.0);
  NumericVector D(n, 0.0);
  
  NumericVector AT(n, 0.0);
  NumericVector DT(n, 0.0);
  
  NumericVector int1(n, 0.0);
  
  for (int i = 0; i < ne; i++) {
    
    Range RN = Range(i*n, (i+1)*n-1);
    NumericVector LaTemp(n);
    computeLAM(LaTemp, lambda[RN], y[RN], n, 0);
    std::copy(LaTemp.begin(), LaTemp.end(), La.begin()+i*n);
    
    for (int j = 0; j < p; j++) {
      Ypre[RN] += X[Range(i*n+j*N, (i+1)*n-1+j*N)] * coef[j];
    }
    
  }
  YpreExp = exp(Ypre);
  
  for (int i = 0; i < ne; i++) {
    
    Range RN = Range(i*n, (i+1)*n-1);
    A += La[RN] * YpreExp[RN];
    for (int j = 0; j < n; j++) {
      if (dn[i*n + j] == 1) {
        B[j] *= lambda[i*n + j] * YpreExp[i*n + j];
      }
    }
    D += dn[RN];
  }
  
  if (frailty == 1 || frailty == 2) {
    
    gsl_set_error_handler_off();
    int status;
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1;
    double result(0.0), error;
    
    if (frailty == 1) {
      F1.function = &logN1int;
    }
    
    if (frailty == 2) {
      F1.function = &InvG1int;
    }
    
    if (frailty == 3) {
      F1.function = &PVF1int;
      kint.mpvf = power;
    }
    
    F1.params = &kint;
    kint.s = tht;
    
    for (int i = 0; i < n; i++) {
      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];
      
      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        
        status = gsl_integration_qags (&F1, 0.001, 10, 0, 1e-7, 1000, w, &result, &error);
        Rcout << "F1: Approximated by Interval\n";
        
      }
      
      int1[i] = result;
      
    }
    
    gsl_integration_workspace_free (w);
    
    return(sum(log(int1)));
    
  }

  if (frailty == 0) {
    AT = 1/tht + A;
    DT = 1/tht + D;

    double l(0.0);
    l -= n*(lgamma(1/tht) + std::log(tht)/tht);
    l += sum(lgamma(DT)) - sum(DT*log(AT));
    l += sum(d * Ypre);

    for (int i = 0; i < N; i++) {
      if (lambda[i] > 0) {
        l += std::log(lambda[i]);
      }
    }

    return l;

  }

  return 0;
}


// [[Rcpp::export]]
double LogLikRE(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0,
                const double& tht0, int frailty, const NumericVector& id, int N, int a, int p, double power) {
  
  
  NumericVector coef = clone(coef0);
  NumericVector lambda = clone(lambda0);
  double tht = tht0;
  
  NumericVector La(N);
  
  NumericVector Ypre(N, 0.0);
  NumericVector YpreExp(N, 0.0);
  intParams kint;
  
  NumericVector A(a, 0.0);
  NumericVector B(a, 1.0);
  NumericVector D(a, 0.0);
  
  NumericVector AT(a, 0.0);
  NumericVector DT(a, 0.0);
  
  NumericVector int1(a, 0.0);
  
  
  int tempID(0);

  computeLAM(La, lambda, y, N, 0);
  
  for (int i = 0; i < p; i++) {
    Ypre += X[Range(i*N, (i+1)*N-1)] * coef[i];
  }
  
  YpreExp = exp(Ypre);
  
  double L0 = 0;
  for (int i = 0; i < N; i++) {
    tempID = id[i];
    if (i == 0 || tempID > id[i-1]) {
      L0 = 0;
    }
    A[tempID] += (La[i] - L0) * YpreExp[i];
    L0 = La[i];
  }
  
  for (int j = 0; j < N; j++) {
    tempID = id[j];
    if (d[j] == 1) {
      B[tempID] *= lambda[j] * YpreExp[j];
      D[tempID]++;
    }
  }
  
  if (frailty == 1 || frailty == 2) {
    
    gsl_set_error_handler_off();
    int status;
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F1;
    double result(0.0), error;
    
    if (frailty == 1) {
      F1.function = &logN1int;
    }
    
    if (frailty == 2) {
      F1.function = &InvG1int;
    }
    
    if (frailty == 3) {
      F1.function = &PVF1int;
      kint.mpvf = power;
    }
    
    F1.params = &kint;
    kint.s = tht;
    
    for (int i = 0; i < a; i++) {
      kint.a = A[i];
      kint.b = B[i];
      kint.d = D[i];
      
      status = gsl_integration_qagiu (&F1, 0, 0, 1e-7, 1000, w, &result, &error);
      
      if (status) {
        
        status = gsl_integration_qags (&F1, 0.001, 10, 0, 1e-7, 1000, w, &result, &error);
        Rcout << "F1: Approximated by Interval\n";
        
      }
      
      int1[i] = result;
      
    }
    
    gsl_integration_workspace_free (w);
    
    return(sum(log(int1)));
  }
  
  if (frailty == 0) {
    AT = 1/tht + A;
    DT = 1/tht + D;
    
    double l(0.0);
    l -= a*(lgamma(1/tht) + std::log(tht)/tht);
    l += sum(lgamma(DT)) - sum(DT*log(AT));
    l += sum(d * Ypre);
    
    for (int i = 0; i < N; i++) {
      if (lambda[i] > 0) {
        l += std::log(lambda[i]);
      }
    }
    
    return l;
    
  } 
  
  return 0.0;
}

