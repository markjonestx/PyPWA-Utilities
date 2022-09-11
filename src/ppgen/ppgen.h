

/* 
 * File:   ppgen.h
 * Author: dennisweygand
 *
 * Created on March 27, 2012, 7:53 PM
 */

#pragma once

#include <cstdlib>
#include <string>
#include <ppgen/ntypes.h>
#include <lib/lorentz.h>
#include <lib/particleType.h>
#include <lib/math/VFour.h>


/* units:  Energy, Mass: GeV
   Angles: rad
   dimensions: cm
*/
#define  BEAM_MASS 0.13957
#define  NEUTRON_MASS 0.93957
#define  TARGET_MASS 0.93827
#define  PROTON_MASS 0.93827
#define  KZERO_MASS 0.49772
#define  KCHARGED_MASS 0.49367
#define  PI_MASS 0.13957

#define  PI0_MASS 0.1349764
#define  OMEGA_MASS 0.78194
#define  ELEC_MASS .000511
#define  GAMMA_MASS 0.0
#define  ETAPRIME_MASS 0.95777

#define  VERTEX_X_CENTER 0.0
#define  VERTEX_Y_CENTER 0.0

#define  VERTEX_X_SIGMA 1.17
#define  VERTEX_Y_SIGMA 1.17


double randm(double, double);

double CMmomentum(double, double, double);

void distribute_vertex(vector3_t *);

float gaussian_gen(void);

double GetBeamP();


void GeneralUsage();

void UsageM1(char *ProcessName);

void UsageM2(char *ProcessName);

void UsageM5(char *ProcessName);

void UsageM6(char *ProcessName);

void UsageM9(char *ProcessName);

void UsageM10(char *ProcessName);

void UsageM11(char *ProcessName);

void UsageM13(char *ProcessName);

void UsageM16(char *ProcessName);

void UsageM26(char *ProcessName);

void Usagecpcmn0N(char *ProcessName);

void Usagecpcm(char *ProcessName);

void Usagec1c2(char *ProcessName);

void MUsage(char *ProcessName);


void pipipi(int argc, char *argv[]);

void pipipi0(int argc, char *argv[]);

void ppipiX_gamma(int, char **);

void npip(int argc, char *argv[]);

void npip_gamma(int argc, char *argv[], int, int *);

void ppi0_gamma(int argc, char *argv[]);

void ppi0(int argc, char *argv[]);

void pippim(int argc, char *argv[]);

void kpkspim(int argc, char *argv[]);

void pipipi0X(int argc, char *argv[]);

void omegapipi(int argc, char *argv[], Particle_t Beam, Particle_t Pi1,
               Particle_t Pi2);

void etaprimepi(int argc, char *argv[], Particle_t Beam, Particle_t Pi);

void omegapipi0(int argc, char *argv[], Particle_t Beam, Particle_t Pi1,
                Particle_t Pi2);

void omegaphi(int argc, char *argv[], Particle_t Beam);

void cpcmn0N(int argc, char *argv[], Particle_t Beam, Particle_t Cplus,
             Particle_t Cminus, Particle_t N0, int decay);

void cpcm(int argc, char *argv[], Particle_t Beam, Particle_t Cplus,
          Particle_t Cminus, Particle_t Recoil);

void cpcmEXP(int argc, char *argv[], Particle_t Beam, Particle_t Cplus,
             Particle_t Cminus, Particle_t Recoil);

void bcpcm(int argc, char *argv[], Particle_t Beam, Particle_t Part1,
           Particle_t Part2, Particle_t Part3);

void
bpn(int argc, char *argv[], Particle_t Beam, Particle_t Part1, Particle_t Part2,
    Particle_t Part3, Particle_t D1,
    Particle_t D2);

void body3(int argc, char *argv[], Particle_t Beam, Particle_t Cplus,
           Particle_t Cminus, Particle_t N0, int decay);

double
Decay2(math::VFour &Rest, double massRest, math::VFour &P1, double massP1, math::VFour &P2,
       double massP2);

void
c1c2(int argc, char **argv, Particle_t Beam, Particle_t Target, Particle_t C1,
     Particle_t P1_1, Particle_t P1_2,
     Particle_t C2, Particle_t P2_1, Particle_t P2_2, Particle_t Recoil);

void
twoBodyDecay(double mass, double mass1, double mass2, math::VFour &p1, math::VFour &p2);

void
doubleDalitz(int argc, char *argv[], Particle_t Beam, Particle_t Pi1,
             Particle_t Pi2, int decay);

void generateBeamMCinsideTarget(vector3_t *, vector3_t *);

double
cosThetaExp(double plab, double beamMass, double targetMass,
            double resonanceMass, double recoilMass, double slope);

float tmin(float plab, float ma, float mb, float mc, float md);

float tmax(float plab, float ma, float mb, float mc, float md);

double
tprimeExp(double plab, double beamMass, double targetMass, double resonanceMass,
          double recoilMass, double slope);

double expDist(double);

double getT(double tMin, double slope);

float pprime(float w, float m, float m3);

float wcm(float p, float m1, float m2);

float e(float p, float m);

float theta(float QSQ, float W, float E);

float eprime(float theta, float Qsq, float E);

double Mass(Particle_t pid);

void PrintParticleID();
