//auxiliary subroutines and functions for describing sets of identical qubits

#pragma once

#include <cmath>

//const double M_PI = 3.1415926535897932385;        //only in Windows

double binomial(const int n, const int k);

inline double pq(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double dpqdtheta1(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double dpqdtheta2(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double dpqdphi1(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double dpqdphi2(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double d2pqdtheta1dtheta1(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double d2pqdtheta1dtheta2(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double d2pqdtheta1dphi1(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double d2pqdtheta1dphi2(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double d2pqdtheta2dtheta1(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double d2pqdtheta2dtheta2(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double d2pqdtheta2dphi1(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
inline double d2pqdtheta2dphi2(double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);

double pt(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double dptdtheta1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double dptdtheta2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double dptdphi1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double dptdphi2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdtheta1dtheta1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdtheta1dtheta2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdtheta1dphi1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdtheta1dphi2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdtheta2dtheta1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdtheta2dtheta2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdtheta2dphi1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdtheta2dphi2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdphi1dtheta1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdphi1dtheta2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdphi1dphi1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdphi1dphi2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdphi2dtheta1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdphi2dtheta2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdphi2dphi1(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);
double d2ptdphi2dphi2(int const nq, int const k, double const p0, double const p1, double const theta1, double const phi1, double const theta2, double const phi2);

inline long double pq(long double const p0, long double const p1, long double const theta1, long double const phi1, long double const theta2, long double const phi2);

void adjust_into_interval(double const a, double const b, double& v);

void ml_estimator_1m(double const p0, double const p1, int const N1, int const nf1, double& thetaf);
void ml_estimator_2m(double const p0, double const p1, int const N1, int const N2, int const nf1, int const nf2, double const thetaf2, double& thetaf, double& phif);
void ml_estimator_3m(double const p0, double const p1, int const N1, int const N2, int const N3, int const nf1, int const nf2, int const nf3,
    double const thetaf2, double const thetaf3, double const phif3, double& thetaf, double& phif);

void direct_inversion_tomography(int const N1, int const N2, int const N3, int const nf1, int const nf2, int const nf3, double& thetaf, double& phif);
