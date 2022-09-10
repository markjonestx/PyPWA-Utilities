#pragma once
#include <iostream>
#include <string>
#include <complex>
#include <cstdio>


typedef enum {
    g_Unknown = 0,
    g_Gamma = 1,
    g_Positron = 2,
    g_Electron = 3,
    g_Pi0 = 7,
    g_PiPlus = 8,
    g_PiMinus = 9,
    g_KLong = 10,
    g_KPlus = 11,
    g_KMinus = 12,
    g_Neutron = 13,
    g_Proton = 14,
    g_AntiProton = 15,
    g_KShort = 16,
    g_Eta = 17,
    g_Lambda = 18,
    g_Deuteron = 45,
    g_omega = 60,
    g_EtaPrime = 61,
    g_phiMeson = 62
} Geant_ID;


inline double tilde(int l) {
    return pow(l + 1, 0.5);
}


#define signof(x) (x<0 ? -1 : 1)
#define MAX(x, y) (x>y ? x : y)
#define MIN(x, y) (x<y ? x : y)


std::complex<double>
D(double alpha, double beta, double gamma, int j, int n, int m);

double clebsch(int j1, int j2, int j3, int m1, int m2, int m3);

double d_jmn_b(int J, int M, int N, double beta);

int fact(int i);

double dfact(double i);

double F(int n, double p);

double lambda(double a, double b, double c);

std::complex<double> q(double M, double m1, double m2);


void addtab();

void subtab();

void ptab();

std::string itos(int);

std::string chargetos(int);

Geant_ID name2id(std::string name, int q);

std::string id2name(Geant_ID);
