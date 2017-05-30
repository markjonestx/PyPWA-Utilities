/* 
 * File:   main.cc
 * Author: dennisweygand
 *
 * Created on March 28, 2012, 7:56 AM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */

/* 
 * File:   ppgen.cc
 * Author: dennisweygand
 *
 * Created on March 27, 2012, 7:47 PM
 */

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <ppgen.h>
#include <event.h>
#include <txtEvent.h>

using namespace std;



extern "C" {
#include <stdio.h>
	   }

using namespace std;


#define SPEED_OF_LIGHT 3.0E8

void pParticle_txt2part(Particle_t type,threeVec v,fourVec p);
void pParticle(Particle_t,threeVec,fourVec);
void pParticleGamp(Particle_t,fourVec);
int  Q(Particle_t  pid);
double brem(double brem0,double brem1);
std::string pid2name(Particle_t  type);

void Usage(char *);

static int runNo = 7000;
static float pbeamZ = 4.0;
static float EbeamZ = 4.0;
static float TLow=0.0,TMax = 0.0;
static float z = 0.0;
int verbose = 0;

static float range0 = 0.0,range1 = 0.0;
static double brem0 = 1.0,brem1 = 4.0;
int UseBrem = 0;
int UseRange = 0;
static float targetZ0 = -9.0,targetZ1 = 9.0;
int UseTarget = 0;
int FlatMass = 0;

int printBeam = 1;
int printAll = 0;
int Isotropic = 0;
int uChannel = 0;
double Slope = -1.0;
int txt2part_style = 0;
int gamp = 0;
int dalitz = 0;

extern particleDataTable PDGtable;



void test100(int argc, char *argv[]);
void UsageM100 (char *ProcessName);
double Qsq(double E,double Ep,double theta);
double Wsq(double E,double Ep,double theta);



int onList(int a,int n,int *l)
{
  int ret = 0;
  if (n) {
    for (int i = 0; i < n; ++i) {
      if (*l++ == a) 
	return(1);
    }
    return(0);
  }
  else
    return(1);
}

int countchar(char *ptr,char c)
{
  int ret = 0;
  while (*ptr) {
    if (*ptr++ == c) 
      ret++;
  }
  return(ret);
}
    

int processList(char *arg,int **alist)
{
  char *word;
  int i;
  int nwords = countchar(arg,',') + 1;
  int *list = new int[nwords];
  word = strtok(arg,"\n,");
  for (i = 0; i < nwords; ++i) {
   list[i] = atoi(word);
   word = strtok(NULL,"\n,");
  }
  *alist = list;
  return(nwords);
}

int main (int argc, char *argv[])
{
  int mode = 0;
  int newargc = 0;
  int nw = 0,*wlist;
  //  ios::sync_with_stdio();
  Particle_t beamid = Gamma;

  char *word;

  char **newargv = (char **) malloc(argc * sizeof(char *));

  
  newargv[newargc++] = argv[0];	   
  
  long int RDMseed = 0;

  //srand48(seed);
  //  srand48((long)time(NULL));

  for (int i = 0; i < argc; ++i) {
    std::cerr << argv[i] << " ";
  }
  std::cerr << std::endl; 

  PDGtable.initialize();

  if (argc <= 1) {
    MUsage (argv[0]);
    exit(1);
  }
  else {
    int usage = 1;
    for (int iarg = 1; iarg < argc; ++iarg) {
      char *ptr = argv[iarg];
      if (*ptr == '-') {
	ptr++;
	switch (*ptr) {
	case 'h':
	  if (usage) {
	    MUsage (argv[0]);
	    std::cerr << "\n\nThe Particle Data Table: " << std::endl;
	    PDGtable.print();
	    exit(0);
	  }
	  else 
	    newargv[newargc++] = argv[iarg];	     
	  break;
	case 'S':
	  RDMseed = atoi(++ptr);
	  break;
	case 'H':
	  if (usage) {
	    MUsage (argv[0]);
	    PrintParticleID();
	    exit(0);
	  }
	  else 
	    newargv[newargc++] = argv[iarg];	     
	  break;

	case 'F':
	  FlatMass = 1;
	  break;

	case 'X':
	  ptr++;
	  word = strtok(ptr,",");
	  TLow = atof(word);
	  word = strtok(NULL," ");
	  TMax = atof(word);
	  std::cerr << "t range: " << TLow << " --> " << TMax << std::endl;
	  break;
	case 'A':
	  // print all particles
	  printAll = 1;
	  break;
	case 'j':
	  // Beam id
	  beamid = (Particle_t) atoi(++ptr);
	  break;
	case 'r':
	  ptr++;
	  word = strtok(ptr,",");
	  range0 = atof(word);
	  word = strtok(NULL," ");
	  range1 = atof(word);
	  std::cerr << "Beam range: " << range0 << " --> " << range1 << std::endl;
	  UseRange = 1;
	  break;
	case 'z':
	  ptr++;
	  if (strlen(ptr)) {
	    word = strtok(ptr,",");
	    targetZ0 = atof(word);
	    word = strtok(NULL," ");
	    targetZ1 = atof(word);
	  }
	  std::cerr << "target range: " << targetZ0 << " --> " << targetZ1 << std::endl;
	  UseTarget = 1;
	  break;
	case 'b':
	  ptr++;
	  word = strtok(ptr,",");
	  brem0 = atof(word);
	  word = strtok(NULL," ");
	  brem1 = atof(word);
	  std::cerr << "Brem range: " << brem0 << " --> " << brem1 << std::endl;
	  UseBrem = 1;
	  break;
	case 'v':
	  verbose = 1;
	  break;
	case 'I':
	  Isotropic = 1;
	  break;
	case 'Z':
	  z = atof(++ptr);
	  break;
	case 'R':
	  runNo = atoi(++ptr);
	  break;
	case 'P':
	  /* beam momentum */
	  pbeamZ = atof(++ptr);
	  break;
	case 'E':
	  /* for gamma beams- primary electron beam energy */
	  EbeamZ = atof(++ptr);
	  break;
	case 't':
	  Slope = atof(++ptr);
	  break;
	case 'M':
	  mode = atoi (++ptr);
	  usage = 0;
	  break;
	  /*	case 'l':
	  nw = processList(++ptr,&wlist);
	  std::cerr << "writing: ";
	  for (int i = 0; i < nw; ++i)
	    std::cerr << wlist[i] << " " ;
	  std::cerr << std::endl;
	  break; */
	case 'T':
	  txt2part_style = 1;
	  break;
	case 'G':
	  gamp = 1;
	  break;
	case 'd':
	  dalitz = 1;
	  break;
	default:
	  newargv[newargc++] = argv[iarg];
	  break;
	}
      }
    }
    if(!RDMseed){
      srand48((long)time(NULL));
    }else{
      srand48(RDMseed);
    }
    
    switch (mode) {
    case 1:
      pipipi (newargc, newargv);
      break;
    case 2:
      pippim(newargc, newargv);
      break;
    case 3:
      npip(newargc,newargv);
      break;
    case 4:
      ppi0(newargc,newargv);
      break;
    case 5:
      kpkspim(newargc,newargv);
      break;
    case 6:
      npip_gamma(newargc,newargv,nw,wlist);
      break;
    case 7:
      cpcmn0N(newargc,newargv,beamid,PiPlus,PiMinus,PiPlus,0);
      break; 
    case 8:
      cpcmn0N(newargc,newargv,beamid,PiPlus,PiMinus,Pi0,1);
      break;
    case 9:
      ppi0_gamma(newargc,newargv);
      break;
    case 10:
      pipipi0(newargc,newargv);
      break;
    case 11:
      //      ppipi_gamma(newargc,newargv);
      cpcm(newargc,newargv,beamid,PiPlus,PiMinus,Proton);
      break;
    case 12:
      ppipiX_gamma(newargc,newargv);
      break;
    case 13:
      pipipi0X(newargc,newargv);
      break;
    case 14:
      //    nKpKm(newargc,newargv);
      cpcm(newargc,newargv,beamid,KPlus,KMinus,Neutron);
      break;
    case 15:
      cpcm(newargc,newargv,beamid,Positron,Electron,Proton);
      break;
    case 16:
      omegapipi(newargc,newargv,beamid,PiPlus,PiMinus);
      break; 
    case 17:
      omegapipi0(newargc,newargv,beamid,PiMinus,Pi0);
      break;
    case 18:
      cpcmn0N(newargc,newargv,beamid,KPlus,KMinus,Pi0,1);
      break;
    case 19:
      cpcmn0N(newargc,newargv,beamid,KPlus,KMinus,Eta,1);
      break;
    case 20:
      cpcm(newargc,newargv,beamid,KMinus,Proton,KPlus);
      break;
     case 21:
      cpcmn0N(newargc,newargv,beamid,PiPlus,PiMinus,Eta,1);
      break; 
    case 22:
      cpcm(newargc,newargv,beamid,PiPlus,PiMinus,Proton);
      break;
    case 23:
      cpcmn0N(newargc,newargv,beamid,PiPlus,PiPlus,PiMinus,0);
      break;
    case 24:
       cpcm(newargc,newargv,beamid,KPlus,KMinus,Proton);
      break;
    case 25:
      cpcm(newargc,newargv,beamid,Proton,AntiProton,Proton);
      break;
    case 26:
      etaprimepi(newargc,newargv,beamid,PiMinus);
      break; 
    case 27:
      cpcmn0N(newargc,newargv,beamid,PiPlus,PiMinus,PiMinus,0);
      break;
    case 28:
      body3(newargc,newargv,beamid,Proton,Proton,AntiProton,0);
      break; 
    case 29:
      cpcmn0N(newargc,newargv,beamid,KPlus,KMinus,PiPlus,0);
      break;
    case 30:
      c1c2(newargc,newargv,beamid,Proton,phiMeson,KPlus,KMinus,Eta,Gamma,Gamma,Proton);
      break;
    case 31:
      bcpcm(newargc,newargv,beamid,KPlus,Proton,KMinus);
      break;  
    case 32:
      bcpcm(newargc,newargv,beamid,PiPlus,Proton,PiMinus);
      break; 
    case 33:
      bcpcm(newargc,newargv,beamid,PiMinus,Proton,PiPlus);
      break; 


    case 34:
       bcpcm(newargc,newargv,beamid,KMinus,Neutron,KPlus);
      break;   

    case 35:
      body3(newargc,newargv,beamid,Proton,PiMinus,PiPlus,0);
      break;

    case 36:
      bpn(newargc,newargv,beamid,KShort,KPlus,Neutron,PiPlus,PiMinus);
      break;

    case 37:
      cpcmEXP(newargc,newargv,beamid,Positron,Electron,Proton);
      break;
    case 38:
      doubleDalitz(newargc,newargv,beamid,Pi0,Pi0,Proton,1);
      break;

    case 39:
      body3(newargc,newargv,beamid,Proton,Pi0,Pi0,0);
      break;

    case 40:
      c1c2(newargc,newargv,beamid,Alpha,Pi0,Gamma,Gamma,Eta,Gamma,Gamma,Alpha);
      break;


   case 41:
      c1c2(newargc,newargv,beamid,Proton,Pi0,Gamma,Gamma,Eta,Gamma,Gamma,Proton);
      break;
   case 42:
       c1c2(newargc,newargv,beamid,Proton,KShort,PiPlus,PiMinus,KShort,PiPlus,PiMinus,Proton);
      break;
   case 43:
       omegaphi(newargc,newargv,beamid);
      break;

    case 100:
      test100(newargc, newargv);
      break;
    case 101:
      bcpcm(newargc,newargv,beamid,KPlus,Proton,KMinus);
      break;

    default:
      MUsage (argv[0]);
      break;
    }
  }
}



// resonance -> pi+ pi-
// isobar1 -> piplus + piminus2

void pippim (int argc, char *argv[])
{
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double masslow = 0.0,
    masshigh = 3.0,
    r2 = 0.0,
    t_min = 0.0,
    t_max,
    slope = 3.0,
    LorentzFactor = 0,
    expt_min = exp (-slope * t_min),
    expt_max,
    lfmax = 0,
    resonance_mass;
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);

  fourVec beam,
    target,
    resonance,
    recoil,
    pip,
    pim;

  vector3_t 
    vbeam,
    pbeam;

  lorentzTransform Boost;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof (ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof (ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'h':
	UsageM2 (argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }
  

  /*
   *-- beam and target in lab frame
   */


  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec (sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) +pow((double) pbeam.z,2.0) + pow (BEAM_MASS, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();
    double pi_pi_threshold = 2 * PI_MASS;

    /* in case it was not set */

    if (masslow == 0.) {
      std::cerr << "pi pi low mass below pi pi threshold- resetting" << std::endl;
      masslow = pi_pi_threshold;
    }
    if(masshigh == 0.)
      masshigh = CMenergy - PROTON_MASS;

    if (masshigh < pi_pi_threshold) {
      std::cerr << "pi pi high mass below pi pi threshold" << std::endl;
      exit(1);
    }

    if (masshigh > CMenergy - PROTON_MASS)
      masshigh = CMenergy - PROTON_MASS;
    if (masslow > CMenergy - PROTON_MASS) {
      std::cerr << "Not enough beam energy... exiting" << std::endl;
      exit(1);
    }

    t_max = pow (CMmomentum (2 * CMenergy, masslow, PROTON_MASS), 2.0);
    expt_max = exp (-slope * t_max);

    /*
     *-- generate the resonance
     */

    if (masshigh < pi_pi_threshold) {
      std::cerr << "pi pi high mass below pi pi_threshold" << std::endl;
      exit (1);
    }

    do {
      resonance_mass = randm (masslow, masshigh);
    } while (resonance_mass < pi_pi_threshold);

    r2 = randm (0.0, 1.0);

    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double expt;
    double t;
    do {
      expt = randm (expt_min, expt_max);
      t = -log (expt) / slope;
      costheta = (CMbeam.t() * resonance_E
		  - 0.5 * (t + pow (PI_MASS, 2.0) + pow (resonance_mass, 2.0)))
	/ (~CMbeam.V() * resonance_p);
    } while (fabs (costheta) > 1.0);
    resonance.polar(resonance_p,acos (costheta), randm (-3.1415, 3.1415));
    resonance.t(resonance_E);

    /*
     *-- recoil particle
     */

    recoil.set(sqrt(resonance.V().lenSq() + pow (PROTON_MASS, 2.0)), zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    double pip_p = CMmomentum (resonance_mass, PI_MASS, PI_MASS);
    pip.polar(pip_p, acos (randm (-0.9999, 0.9999)),randm (-M_PI,M_PI));
    pip.t(sqrt (pip.V().lenSq () + pow (PI_MASS, 2.0)));

    pim.set(sqrt(pip.V().lenSq() +  pow (PI_MASS, 2.0)),zeroVec - pip.V());

    /*
     *  compute lorentz factor
     */
    LorentzFactor = resonance_p * (pip_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	//	fourVec tmpvec;
	/*
	 *  transform all 4-vectors back to lab frame
	 */
	fourVec tmp;

	// boost from pi_pi rest frame to CM
	tmp.V() = zeroVec - resonance.V();
	tmp.t(resonance.t());
	Boost.set (tmp);
	pip *= Boost;
	pim *= Boost;

	// boost from CM to target rest frame (lab)
	Boost.set (CMtarget);
	resonance *= Boost;
	recoil *= Boost;
	pip *= Boost;
	pim *= Boost;

	/*
	 *  generate vertices
	 */
	// e852 values
	threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	//	threeVec decay;

	if (Print) {
	  // std::cerr << "\n\n*** New Event\n";
	  //	  std::cerr << "Beam:\n  " << beam;
	  // std::cerr << "Beam in CM:\n  " << CMbeam;
	  // std::cerr << "Target in CM:\n  " << CMtarget;
	  // std::cerr << "Resonance\n";
	  // std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  // std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  // std::cerr << "  expt: " << expt << "\n";
	  // std::cerr << "  t: " << t << "\n";
	  // std::cerr << "  Resonance:\n     " << resonance;
	  // std::cerr << "recoil: \n  " << recoil;
	  // std::cerr << "pip\n";
	  // std::cerr << "  pip mass: " << PI_MASS << "\n";
	  // std::cerr << "  pip momentum: " << pip_p << "\n";
	  // std::cerr << "  pip:\n    " << pip;
	  // std::cerr << "pim :\n  " << pim;
	  // std::cerr << "vertices:\n";
	  // std::cerr << "  prod: " << production;
	  // std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cerr << pow (~(pip + pim), 2.0) << std::endl;
	}
	/*
	 *  write event
	 */
	/*	  pParticle(Electron,production,beam); */
	pParticle(PiPlus,production,pip);
	pParticle(PiMinus,production,pim);
	/*	  pParticle(Proton,production,recoil); */
	maxevents--;
      }
    }
  }
}

/*------------- end of pi+ pi- -------------------------------*/

void nKpKm(int argc, char *argv[])
{
  float beamMass = PI_MASS;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
    double masslow = 2 * KCHARGED_MASS,
    masshigh = 3.0,
    t_min = 0.0,
    t_max,
    slope = 3.0,
    LorentzFactor = 0,
    expt_min = exp (-slope * t_min),
    expt_max,
    lfmax = 0,
    resonance_mass;
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);

  fourVec beam,
    target,
    resonance,
    recoil,
    Kp,
    Km;

  vector3_t 
    vbeam,
    pbeam;

  int printBaryon = 0;

  lorentzTransform Boost;

    for (int iarg = 1; iarg < argc; ++iarg) {
      char *ptr = argv[iarg];
      if (*ptr == '-') {
	ptr++;
	switch (*ptr) {
	case 'B':
	  printBaryon = 1;
	  break;
	case 'm':
	  ptr++;
	  maxevents = atoi (ptr);
	  std::cerr << "maxevents: " << maxevents << "\n";
	  break;
	case 'l':
	  ptr++;
	  lfevents = atoi (ptr);
	  std::cerr << "lfevents: " << lfevents << "\n";
	  break;
	case 'L':
	  ptr++;
	  masslow = atof (ptr);
	  break;
	case 'U':
	  ptr++;
	  masshigh = atof (ptr);
	  break;
	case 'p':
	  Print = 1;
	  break;
	case 'M':
	  break;
	case 'b':
	  beamMass = atof(++ptr);
	  break;
	case 'h':
	  UsageM14 (argv[0]);
	  return;
	default:
	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
	  break;
	}
      }
    }
  

  /*
   *-- beam and target in lab frame
   */

    slope = (Slope > 0.0) ? Slope : slope; 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec (sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) +pow((double) pbeam.z,2.0) + pow((double)beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();
   double K_K_threshold = 2 * KCHARGED_MASS;

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PROTON_MASS;

    if (masshigh < K_K_threshold) {
      std::cerr << "K K high mass below K K threshold" << std::endl;
      exit(1);
    }

    if (masshigh > CMenergy - PROTON_MASS)
      masshigh = CMenergy - PROTON_MASS;
    if (masslow > CMenergy - PROTON_MASS) {
      std::cerr << "Not enough beam energy... exiting" << std::endl;
      exit(1);
    }

    t_max = pow (CMmomentum (2 * CMenergy, masslow, PROTON_MASS), 2.0);
    expt_max = exp (-slope * t_max);

    /*
     *-- generate the resonance
     */

    if (masshigh < K_K_threshold) {
      std::cerr << "K K high mass below K K_threshold" << std::endl;
      exit (1);
    }

    do {
      resonance_mass = randm (masslow, masshigh);
    } while (resonance_mass < K_K_threshold);


    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double expt;
    double t;
    do {
      expt = randm (expt_min, expt_max);
      t = -log (expt) / slope;
      costheta = (CMbeam.t() * resonance_E
		  - 0.5 * (t + pow((double) beamMass, 2.0) + pow (resonance_mass, 2.0)))
	/ (~CMbeam.V() * resonance_p);
    } while (fabs (costheta) > 1.0);
    resonance.polar(resonance_p,acos (costheta), randm (-M_PI, M_PI));
    resonance.t(resonance_E);

    /*
     *-- recoil particle
     */

    recoil.set(sqrt(resonance.V().lenSq() + pow (PROTON_MASS, 2.0)), zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    double Kp_p = CMmomentum (resonance_mass, KCHARGED_MASS, KCHARGED_MASS);
    Kp.polar(Kp_p, acos (randm (-1,1)),randm (-M_PI,M_PI));
    Kp.t(sqrt (Kp.V().lenSq () + pow (KCHARGED_MASS, 2.0)));

    Km.set(sqrt(Kp.V().lenSq() +  pow (KCHARGED_MASS, 2.0)),zeroVec - Kp.V());

    /*
     *  compute lorentz factor
     */
    LorentzFactor = resonance_p * (Kp_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	//	fourVec tmpvec;
	/*
	 *  transform all 4-vectors back to lab frame
	 */
	fourVec tmp;

	// boost from K_K rest frame to CM
	tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set (tmp);
	Kp *= Boost;
	Km *= Boost;

	// boost from CM to target rest frame (lab)
	Boost.set (CMtarget);
	resonance *= Boost;
	recoil *= Boost;
	Kp *= Boost;
	Km *= Boost;

	/*
	 *  generate vertices
	 */
	// e852 values
	threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	//	threeVec decay;

	if (Print) {
	  // std::cerr << "\n\n*** New Event\n";
	  //	  std::cerr << "Beam:\n  " << beam;
	  // std::cerr << "Beam in CM:\n  " << CMbeam;
	  // std::cerr << "Target in CM:\n  " << CMtarget;
	  // std::cerr << "Resonance\n";
	  // std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  // std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  // std::cerr << "  expt: " << expt << "\n";
	  // std::cerr << "  t: " << t << "\n";
	  // std::cerr << "  Resonance:\n     " << resonance;
	  // std::cerr << "recoil: \n  " << recoil;
	  // std::cerr << "Kp\n";
	  // std::cerr << "  Kp mass: " << KCHARGED_MASS << "\n";
	  // std::cerr << "  Kp momentum: " << Kp_p << "\n";
	  // std::cerr << "  Kp:\n    " << Kp;
	  // std::cerr << "Km :\n  " << Km;
	  // std::cerr << "vertices:\n";
	  // std::cerr << "  prod: " << production;
	  // std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cerr << pow (~(Kp + Km), 2.0) << std::endl;
	}
	/*
	 *  write event
	 */
	if (txt2part_style) {
	  std::cout << "3" << std::endl;
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;
	  pParticle_txt2part(KPlus,production,Kp);
	  pParticle_txt2part(KMinus,production,Km);
	  pParticle_txt2part(Neutron,production,recoil);
	}
	else {
	  pParticle(KPlus,production,Kp);
	  pParticle(KMinus,production,Km);
	  if (printBaryon)
	    pParticle(Neutron,production,recoil);
	}

	maxevents--;
      }
    }
  }
}

// resonance -> isobar1 + piminus1
// isobar1 -> piplus + piminus2

void pipipi(int argc, char *argv[])
{
  int icount = 0;
  int Print = 0,
  maxevents = 9999999,
  nevents = 0,
    lfevents = 5000;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    slope = 10.0,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    piminus1,
    piplus,
    piminus2,
    isobar1;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = BEAM_MASS;
  int printBaryon = 0;

  float tMin;

    for (int iarg = 1; iarg < argc; ++iarg) {
      char *ptr = argv[iarg];
      if (*ptr == '-') {
	ptr++;
	switch (*ptr) {
	case 'm':
	  ptr++;
	  maxevents = atoi (ptr);
	  std::cerr << "maxevents: " << maxevents << "\n";
	  break;
	case 'l':
	  ptr++;
	  lfevents = atoi (ptr);
	  std::cerr << "lfevents: " << lfevents << "\n";
	  break;
	case 'L':
	  ptr++;
	  masslow = atof (ptr);
	  break;
	case 'U':
	  ptr++;
	  masshigh = atof (ptr);
	  break;
	case 'p':
	  Print = 1;
	  break;
	case 'M':
	  break;
	case 'b':
	  beamMass = atof(++ptr);
	  break;
	case 'B':
	  printBaryon = 1;
	  break;
	case 'h':
	  UsageM1 (argv[0]);
	  return;
	default
:	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
	  break;
	}
      }
    }


  double pi_pi_pi_threshold = 3 * PI_MASS;
  double pi_pi_threshold = 2 * PI_MASS;
  if(masslow < pi_pi_pi_threshold)
    masslow = pi_pi_pi_threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PROTON_MASS;
    if (masshigh < pi_pi_pi_threshold) {
      std::cerr << "pi_pi_pi high mass below pi_pi_pi_threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < pi_pi_pi_threshold);
    isobar1_mass = randm(pi_pi_threshold, (resonance_mass - PI_MASS));

    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,NEUTRON_MASS);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;
    icount = 0;
      do {
	if((icount++) > 100000) {
	  std::cerr << "exiting... too many counts" << std::endl;
	  exit(0);
	}

      t = randm(0., t_max);
      //     expt_max = t_max;
    } while( randm(0.,expt_max) > t*exp(-slope*t) );
    costheta = 1. - 2.*t/t_max;
    //   costheta = 1.0;

    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, PI_MASS);
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

     isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    piminus1.set(sqrt(isobar1.V().lenSq() + pow (PI_MASS, 2.0)),zeroVec -isobar1.V());




    /*
     *  now do decay in isobar1 rest frame
     */
  // pi pi case
    double piplus_p = CMmomentum (isobar1_mass, PI_MASS,PI_MASS);

    piplus.polar( piplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   piplus.polar( piplus_p,acos (1.0),0.0);
    piplus.t(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)));
    piminus2.set(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)),zeroVec - piplus.V());


    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar1_p * piplus_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
         /* transform all 4-vectors back to lab frame */
        fourVec tmp;

        // boost from pi_pi(rho) rest frame to pi_pi_pi(a2/M) rest frame

        tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        piplus = Boost * piplus;
        piminus2 = Boost * piminus2;

        // boost from pi_pi_pi rest frame to CM
        tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);

        isobar1 = Boost * isobar1;
        piplus = Boost * piplus;
        piminus2 = Boost * piminus2;
        piminus1 = Boost * piminus1;

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        isobar1 = Boost * isobar1;
        piplus = Boost * piplus;
        piminus2 = Boost * piminus2;
        piminus1 = Boost * piminus1;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

	if (Print) {
	   std::cerr << "\n\n*** New Event\n";
	   std::cerr << "Beam:\n  ";
	   beam.print();
	   std::cerr << "Beam in CM:\n  ";
	   CMbeam.print();
	   std::cerr << "Target in CM:\n  " ;
	   CMtarget.print();
	   std::cerr << "Resonance\n";
	   std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	   std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	   std::cerr << "  t: " << t << "\n";
	   std::cerr << "  Resonance:\n "    ;
	   resonance.print();
	   std::cerr << "recoil: \n  ";
	   recoil.print();
	   std::cerr << "isobar1\n";
	   std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	   std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	   std::cerr << "  isobar1:\n    ";
	   isobar1.print();
	   std::cerr << "pi+ :\n  ";
	   piplus.print();
	   std::cerr << "pi-2:\n  ";
	   piminus2.print();
	   std::cerr << "pi-1:\n  ";
	   piminus1.print();
	   //	   std::cerr << "vertices:\n";
	   //	   std::cerr << "  prod: " << production;
	   std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	   std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	  if (dalitz) {
	    std::cerr << pow (~(piminus2 + piplus), 2.0) << " "
	      << pow (~(piminus2 + piminus1), 2.0) << " "
		<< pow (~(piplus + piminus1), 2.0) << std::endl;
	  }
	/*
	 *  write event
	 */

	  pParticle(PiMinus,production,piplus);
	  pParticle(PiPlus,production,piminus1);
	  pParticle(PiPlus,production,piminus2);
	  if (printBaryon)
	    pParticle(Proton,production,recoil);
	  maxevents--;
      }
    }
  }
}

/*----------------End of 3 pi-----------------------------------*/



void pipipi0 (int argc, char *argv[])
{
  int icount = 0;
  int printGamma = 0;
  int Print = 0,
  maxevents = 9999999,
  nevents = 0,
  lfevents = 5000;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    slope = 10.0,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    elec,
    resonance,
    recoil,
    piminus,
    piplus,
    pizero,
    isobar1;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = BEAM_MASS;
  int printBaryon = 0;

  float tMin;

    for (int iarg = 1; iarg < argc; ++iarg) {
      char *ptr = argv[iarg];
      if (*ptr == '-') {
	ptr++;
	switch (*ptr) {
	case 'm':
	  ptr++;
	  maxevents = atoi (ptr);
	  std::cerr << "maxevents: " << maxevents << "\n";
	  break;
	case 'l':
	  ptr++;
	  lfevents = atoi (ptr);
	  std::cerr << "lfevents: " << lfevents << "\n";
	  break;
	case 'L':
	  ptr++;
	  masslow = atof (ptr);
	  break;
	case 'U':
	  ptr++;
	  masshigh = atof (ptr);
	  break;
	case 'p':
	  Print = 1;
	  break;
	case 'M':
	  break;
	case 'b':
	  beamMass = atof(++ptr);
	  break;
	case 'B':
	  printBaryon = 1;
	  break;
	case 'h':
	  UsageM10(argv[0]);
	  return;
	default:
	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
	  break;
	}
      }
    }


  double pi_pi_pi_threshold = 3 * PI_MASS;
  double pi_pi_threshold = 2 * PI_MASS;
  if(masslow < pi_pi_pi_threshold)
    masslow = pi_pi_pi_threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PROTON_MASS;
    if (masshigh < pi_pi_pi_threshold) {
      std::cerr << "pi_pi_pi high mass below pi_pi_pi_threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < pi_pi_pi_threshold);
    isobar1_mass = randm(pi_pi_threshold, (resonance_mass - PI_MASS));

    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,NEUTRON_MASS);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;
    icount = 0;
      do {
	if((icount++) > 100000) {
	  std::cerr << "exiting... too many counts" << std::endl;
	  exit(0);
	}

      t = randm(0., t_max);
      //     expt_max = t_max;
    } while( randm(0.,expt_max) > t*exp(-slope*t) );
    costheta = 1. - 2.*t/t_max;
    //   costheta = 1.0;

    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, PI_MASS);
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

     isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    piminus.set(sqrt(isobar1.V().lenSq() + pow (PI_MASS, 2.0)),zeroVec -isobar1.V());




    /*
     *  now do decay in isobar1 rest frame
     */
  // pi pi case
    double piplus_p = CMmomentum (isobar1_mass, PI_MASS,PI_MASS);

    piplus.polar( piplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   piplus.polar( piplus_p,acos (1.0),0.0);
    piplus.t(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)));
    pizero.set(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)),zeroVec - piplus.V());


    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar1_p * piplus_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
         /* transform all 4-vectors back to lab frame */
        fourVec tmp;

        // boost from pi_pi(rho) rest frame to pi_pi_pi(a2/M) rest frame

        tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        piplus = Boost * piplus;
        pizero = Boost * pizero;

        // boost from pi_pi_pi rest frame to CM
        tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);

        isobar1 = Boost * isobar1;
        piplus = Boost * piplus;
        pizero = Boost * pizero;
        piminus = Boost * piminus;

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        isobar1 = Boost * isobar1;
        piplus = Boost * piplus;
        pizero = Boost * pizero;
        piminus = Boost * piminus;


	elec = beam + target - recoil - piminus - piplus - pizero;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

	if (Print) {
	   std::cerr << "\n\n*** New Event\n";
	   std::cerr << "Beam:\n  ";
	   beam.print();
	   std::cerr << "Beam in CM:\n  ";
	   CMbeam.print();
	   std::cerr << "Target in CM:\n  " ;
	   CMtarget.print();
	   std::cerr << "Resonance\n";
	   std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	   std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	   std::cerr << "  t: " << t << "\n";
	   std::cerr << "  Resonance:\n "    ;
	   resonance.print();
	   std::cerr << "recoil: \n  ";
	   recoil.print();
	   std::cerr << "isobar1\n";
	   std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	   std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	   std::cerr << "  isobar1:\n    ";
	   isobar1.print();
	   std::cerr << "pi+ :\n  ";
	   piplus.print();
	   std::cerr << "pi-2:\n  ";
	   pizero.print();
	   std::cerr << "pi-1:\n  ";
	   piminus.print();
	   //	   std::cerr << "vertices:\n";
	   //	   std::cerr << "  prod: " << production;
	   std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	   std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	  if (dalitz) {
	    std::cerr << pow (~(pizero + piplus), 2.0) << " "
	      << pow (~(pizero + piminus), 2.0) << " "
		<< pow (~(piplus + piminus), 2.0) << std::endl;
	  }
	  else if (txt2part_style) {
	    if (printGamma) {
	      std::cout << "6" << std::endl;
	    } 
	    else { 
	      std::cout << "5" << std::endl;
	    }
	    std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;
	    pParticle_txt2part(Electron,production,elec);
	    pParticle_txt2part(Proton, production, recoil);
	    pParticle_txt2part(PiMinus, production, piminus);
	    pParticle_txt2part(PiPlus, production, piplus);
	    if (printGamma){
	      //      pParticle_txt2part(Gamma, production, gamma1);
	      //  pParticle_txt2part(Gamma, production, gamma2);
	    }
	    else{
	      pParticle_txt2part(Pi0,production,pizero);
	    }
	  }
	  else {
	/*
	 *  write event
	 */
	  pParticle(Electron,production,elec);
	  pParticle(PiMinus,production,piplus);
	  pParticle(PiPlus,production,piminus);
	  pParticle(Pi0,production,pizero);
	  if (printBaryon)
	    pParticle(Proton,production,recoil);
	  }
	  maxevents--;
      }
    }
  }
}
void pipipi0X(int argc, char *argv[])
{
  int icount = 0;
  int printGamma = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    slope = 10.0,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    virtualBeam,
    W,
    electron,
    resonance,
    recoil,
    piminus,
    piplus,
    pizero,
    isobar1;
  lorentzTransform Boost;

  fourVec CMbeam;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = BEAM_MASS;
  int printBaryon = 0;

  double qlo = 0,qhi = 4,wlo = 3 * PI_MASS + PROTON_MASS,whi = 4,q,thet,w;
  double ep;
  double px,py,pz,phi;

  float tMin;
  float setMassHigh = 0.0;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof (ptr);
	break;
      case 'U':
	ptr++;
	setMassHigh = atof (ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'q':
	qlo  = atof(++ptr);
	break;
      case 'Q':
	qhi  = atof(++ptr);
	break;
      case 'w':
	wlo  = atof(++ptr);
	break;
      case 'W':
	whi  = atof(++ptr);
	break;

      case 'h':
	UsageM13(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double pi_pi_pi_threshold = 3 * PI_MASS;
  double pi_pi_threshold = 2 * PI_MASS;
  if(masslow < pi_pi_pi_threshold)
    masslow = pi_pi_pi_threshold;

 

  while (maxevents) {

    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);


    /*
     *-- beam and target in lab frame
     */

    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /* get a w and a qsq */

    do {
      w = randm(wlo,whi);
      q = randm(qlo,qhi);
      thet = theta(q,w,beam.t());
    } while (thet < -500.0);
    ep = eprime(thet,q,beam.t());
    phi = randm(0,2.0 * M_PI);
    pz = cos(thet) * ep;
    px = sin(thet) * cos(phi) * ep;
    py = sin(thet) * sin(phi) * ep;
      
    electron = fourVec( sqrt(ep * ep + pow((double) beamMass, 2.0)),threeVec(px,py,pz));

    // the W system is left- W in lab

    W = beam + target  - electron;
 
    double CMenergy = W.t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */

    masshigh = (setMassHigh < pi_pi_pi_threshold) ? CMenergy - PROTON_MASS : setMassHigh;

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < pi_pi_pi_threshold);
    isobar1_mass = randm(pi_pi_threshold, (resonance_mass - PI_MASS));

    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    // Only isotropic right now
    costheta = randm(-1.0,1.0);


    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, PI_MASS);
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

    isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    piminus.set(sqrt(isobar1.V().lenSq() + pow (PI_MASS, 2.0)),zeroVec -isobar1.V());




    /*
     *  now do decay in isobar1 rest frame
     */
    // pi pi case
    double piplus_p = CMmomentum (isobar1_mass, PI_MASS,PI_MASS);

    piplus.polar( piplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   piplus.polar( piplus_p,acos (1.0),0.0);
    piplus.t(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)));
    pizero.set(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)),zeroVec - piplus.V());


    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar1_p * piplus_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

        // boost from pi_pi(rho) rest frame to pi_pi_pi(a2/M) rest frame

        tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        piplus = Boost * piplus;
        pizero = Boost * pizero;

        // boost to W 
        tmp.set(W.t(),zeroVec - W.V());
	Boost.set(tmp);

        isobar1 = Boost * isobar1;
        piplus = Boost * piplus;
        pizero = Boost * pizero;
        piminus = Boost * piminus;

	// boost from W to lab

	Boost.set(target);
	W *= Boost;
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        isobar1 = Boost * isobar1;
        piplus = Boost * piplus;
        pizero = Boost * pizero;
        piminus = Boost * piminus;

  
        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "W:\n  ";
	  W.print();
	  std::cerr << "Beam  in CM:\n  " ;
	  CMbeam.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "isobar1\n";
	  std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	  std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	  std::cerr << "  isobar1:\n    ";
	  isobar1.print();
	  std::cerr << "pi+ :\n  ";
	  piplus.print();
	  std::cerr << "pi-2:\n  ";
	  pizero.print();
	  std::cerr << "pi-1:\n  ";
	  piminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cerr << pow (~(pizero + piplus), 2.0) << " "
	       << pow (~(pizero + piminus), 2.0) << " "
	       << pow (~(piplus + piminus), 2.0) << std::endl;
	}
	else if (txt2part_style) {
	  if (printGamma) {
	    std::cout << "6" << std::endl;
	  } 
	  else { 
	    std::cout << "5" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;
	  pParticle_txt2part(Electron,production,electron);
	  pParticle_txt2part(Proton, production, recoil);
	  pParticle_txt2part(PiMinus, production, piminus);
	  pParticle_txt2part(PiPlus, production, piplus);
	  if (printGamma){
	    //      pParticle_txt2part(Gamma, production, gamma1);
	    //  pParticle_txt2part(Gamma, production, gamma2);
	  }
	  else{
	    pParticle_txt2part(Pi0,production,pizero);
	  }
	}
	else {
	  /*
	   *  write event
	   */
	  pParticle(Electron,production,electron);
	  pParticle(PiMinus,production,piplus);
	  pParticle(PiPlus,production,piminus);
	  pParticle(Pi0,production,pizero);
	  if (printBaryon)
	    pParticle(Proton,production,recoil);
	}
	maxevents--;
      }
    }
  }
}

/*----------------End of 3 pi-----------------------------------*/



// resonance -> pi+ pi-
// isobar1 -> piplus + piminus2

void npip (int argc, char *argv[])
{
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double masslow = 2 * PI_MASS,
    masshigh = 3.0,
    r2 = 0.0,
    t_min = 0.0,
    t_max,
    slope = 3.0,
    LorentzFactor = 0,
    expt_min = exp (-slope * t_min),
    expt_max,
    lfmax = 0,
    resonance_mass;
  fourVec beam,
    target,
    electron,
    resonance,
    neutron,
    pip;

  vector3_t 
    vbeam,
    pbeam;

  lorentzTransform Boost;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof (ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof (ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;

      case 'h':
	UsageM2 (argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }
  

  /*
   *-- beam and target in lab frame
   */


  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec (sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) +
                                                       pow((double) pbeam.z,2.0) + pow (BEAM_MASS, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));

    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));
    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();
    double n_pi_threshold = NEUTRON_MASS + PI_MASS;

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PI_MASS;

    if (masshigh < n_pi_threshold) {
      std::cerr << "n pi high mass below pi pi threshold" << std::endl;
      exit(1);
    }
    if (masshigh > CMenergy - PI_MASS)
      masshigh = CMenergy - PI_MASS;
    if (masslow > CMenergy - PI_MASS) {
      std::cerr << "Not enough beam energy... exiting" << std::endl;
      exit(1);
    }

    t_max = pow (CMmomentum (2 * CMenergy, masslow, PROTON_MASS), 2.0);
    expt_max = exp (-slope * t_max);

    /*
     *-- generate the resonance
     */
    do {
      resonance_mass = randm (masslow, masshigh);
    } while (resonance_mass < n_pi_threshold);

    r2 = randm (0.0, 1.0);

    double resonance_p = CMmomentum (CMenergy, resonance_mass, PI_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double expt;
    double t;
    double tmin; 
    tmin = pow(resonance_mass,2.0) + pow(PROTON_MASS,2) - 2.0 * CMtarget.t() * resonance_E + 2.0 * ~CMtarget.V() * resonance_p;
    expt_min = exp (-slope * tmin);

    do {
      expt = randm (expt_min, expt_max);
      t = -log (expt) / slope;
      /* t -s really -t */
      costheta =- (CMtarget.t() * resonance_E
		   - 0.5 * (t + pow (PROTON_MASS, 2.0) + pow (resonance_mass, 2.0)))
	/ (~CMtarget.V() * resonance_p);
    } while (fabs (costheta) > 1.0);

    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);



    /*
     *  now do decay in resonance rest frame
     */
    double pip_p = CMmomentum (resonance_mass, NEUTRON_MASS, PI_MASS);


    pip.polar( pip_p, acos (randm (-0.9999, 0.9999)),randm(-M_PI,M_PI));
    pip.t(sqrt (pip.V().lenSq () + pow (PI_MASS, 2.0)));
    neutron.V() = threeVec (0.0, 0.0, 0.0) - pip.V();
    neutron.t(sqrt (neutron.V().lenSq () + pow (NEUTRON_MASS, 2.0)));

    /*
     *  compute lorentz factor
     */
    LorentzFactor = resonance_p * (pip_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	//	fourVec tmpvec;
	/*
	 *  transform all 4-vectors back to lab frame
	 */
	fourVec tmp;
	threeVec zeroVec = threeVec (0.0, 0.0, 0.0);

	// boost from n_pi rest frame to CM
	tmp.V() = zeroVec - resonance.V();
	tmp.t(resonance.t());
	Boost.set (tmp);
	pip *= Boost;
	neutron *= Boost;

	// boost from CM to target rest frame (lab)
	Boost.set (CMtarget);
	resonance *= Boost;
	pip *= Boost;
	neutron *= Boost;

	/*
	 *  generate vertices
	 */
	// e852 values
	threeVec production = threeVec (0.0, 0.0, 0.0);
	//	threeVec decay;
	electron = beam + target -neutron - pip;
	if (Print) {
	  // std::cerr << "\n\n*** New Event\n";
	  // std::cerr << "Beam:\n  " << beam;
	  // std::cerr << "Beam in CM:\n  " << CMbeam;
	  // std::cerr << "Target in CM:\n  " << CMtarget;
	  // std::cerr << "Resonance\n";
	  // std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  // std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  // std::cerr << "  expt: " << expt << "\n";
	  // std::cerr << "  t: " << t << "\n";
	  // std::cerr << " Qsq: "  << -(beam - electron).lenSq() << std::endl; 
	  // std::cerr << " t: " << (target - neutron - pip).lenSq() << std::endl;
	  // std::cerr << "  Resonance:\n     " << resonance;
	  // std::cerr << "pip\n";
	  // std::cerr << "  pip mass: " << PI_MASS << "\n";
	  // std::cerr << "  pip momentum: " << pip_p << "\n";
	  // std::cerr << "  pip:\n    " << pip;
	  // std::cerr << "  neutron :\n  " << neutron;
	  // std::cerr << "  electron :\n  " << electron;
	  // std::cerr << "vertices:\n";
	  // std::cerr << "  prod: " << production;
	  // std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cerr << pow (~(pip + neutron), 2.0) << std::endl;
	}
	/*
	 *  write event
	 */
	pParticle(Electron,production,electron);
	pParticle(PiPlus,production,pip);
	pParticle(Neutron,production,neutron);
	maxevents--;
      }
    }
  }
}



void npip_gamma (int argc, char *argv[],int nw,int *wlist)
{
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
    double mass,
    resonance_mass;
  double beamMass = 0.0;
  fourVec beam,
    target,
    resonance,
    neutron,
    pip;

  vector3_t 
    vbeam,
    pbeam;

  lorentzTransform Boost;

    for (int iarg = 1; iarg < argc; ++iarg) {
      char *ptr = argv[iarg];
      if (*ptr == '-') {
	ptr++;
	switch (*ptr) {
	case 'm':
	  ptr++;
	  maxevents = atoi (ptr);
	  std::cerr << "maxevents: " << maxevents << "\n";
	  break;
	case 'l':
	  ptr++;
	  lfevents = atoi (ptr);
	  std::cerr << "lfevents: " << lfevents << "\n";
	  break;
	case 'p':
	  Print = 1;
	  break;
	case 'M':
	  break;

	case 'h':
	  UsageM6 (argv[0]);
	  return;
	default:
	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
	  break;
	}
      }
    }
  

  /*
   *-- beam and target in lab frame
   */


  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec (sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) +
                                                       pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
     target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));
    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();
   double n_pi_threshold = NEUTRON_MASS + PI_MASS;
   threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
      mass = CMenergy;


    /*
     *-- generate the resonance
     */

      resonance_mass = mass;

      if (resonance_mass < NEUTRON_MASS + PI_MASS) {
	std::cerr << "Not enough CM energy...exiting." << std::endl;
	exit(1);
      }


    double resonance_p = 0.0;
    double resonance_E = resonance_mass;
    double costheta;
    double expt;
    double t;

    resonance.set(resonance_mass,zeroVec);

    /*
     *  now do decay in resonance rest frame
     */
    double pip_p = CMmomentum (resonance_mass, NEUTRON_MASS, PI_MASS);

    pip.polar( pip_p, acos (randm (-0.9999, 0.9999)),randm(-M_PI,M_PI));
    pip.t(sqrt (pip.V().lenSq () + pow (PI_MASS, 2.0)));
    neutron.set(sqrt (pip.V().lenSq() + pow (NEUTRON_MASS, 2.0)),zeroVec - pip.V());
 
  

	/*
	 *  transform all 4-vectors back to lab frame
	 */
	fourVec tmp;


	// boost from CM to target rest frame (lab)
	Boost.set (CMtarget);
	resonance *= Boost;
	pip *= Boost;
	neutron *= Boost;

	/*
	 *  generate vertices
	 */
	threeVec production = threeVec (0.0, 0.0, 0.0);
	//	threeVec decay;
	if (Print) {
	  // std::cerr << "\n\n*** New Event\n";
	  // std::cerr << "Beam:\n  " << beam;
	  // std::cerr << "Beam in CM:\n  " << CMbeam;
	  // std::cerr << "Target in CM:\n  " << CMtarget;
	  // std::cerr << "Resonance\n";
	  // std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  // std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  // std::cerr << "  Resonance:\n     " << resonance;
	  // std::cerr << "  pip\n";
	  // std::cerr << "  pip mass: " << PI_MASS << "\n";
	  // std::cerr << "  pip momentum: " << pip_p << "\n";
	  // std::cerr << "  pip:\n    " << pip;
	  // std::cerr << "  neutron :\n  " << neutron;
	  // std::cerr << "vertices:\n";
	  // std::cerr << "  prod: " << production;
	}

	/*
	 *  write event
	 */

	if (txt2part_style) {
	  std::cout << "2" << std::endl;
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;
	    pParticle_txt2part(Neutron, production, neutron);
	    pParticle_txt2part(PiPlus, production, pip);
	}

	else {

	if (onList((int) PiPlus,nw,wlist))
	    pParticle(PiPlus,production,pip);
	if (onList((int) Neutron,nw,wlist))
	  pParticle(Neutron,production,neutron);
	}
	maxevents--;
  }
}










void ppi0_gamma (int argc, char *argv[])
{
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
    int debug = 0;
  double mass,
    resonance_mass;
  fourVec beam,
    target,
    resonance,
    proton,
    gamma1,
    gamma2,
    pi0;

  vector3_t 
    vbeam,
    pbeam;
 char *word;

  int printBaryon = 0;
  int printGamma = 0;
  double beamMass = 0.0;

  lorentzTransform Boost;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'p':
	Print = 1;
	break;
      case 'D':
	debug = 1;
	break;
      case 'M':
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'g':
	printGamma = 1;
	break;
      case 'h':
	UsageM9 (argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }
  

  /*
   *-- beam and target in lab frame
   */


  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec (sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) +
			 pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));
    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();
    double p_pi0_threshold = PROTON_MASS + PI0_MASS;
    threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
    mass = CMenergy;


    /*
     *-- generate the resonance
     */

    resonance_mass = mass;

    if (resonance_mass < PROTON_MASS + PI0_MASS) {
      std::cerr << "Not enough CM energy...exiting." << std::endl;
      exit(1);
    }


    double resonance_p = 0.0;
    double resonance_E = resonance_mass;
    double costheta;
    double expt;
    double t;

    resonance.set(resonance_mass,zeroVec);

    /*
     *  now do decay in resonance rest frame
     */
    double pi0_p = CMmomentum (resonance_mass, PROTON_MASS, PI0_MASS);

    pi0.polar( pi0_p, acos (randm (-0.9999, 0.9999)),randm(-M_PI,M_PI));
    pi0.t(sqrt (pi0.V().lenSq () + pow (PI0_MASS, 2.0)));
    proton.set(sqrt (pi0.V().lenSq() + pow (PROTON_MASS, 2.0)),zeroVec - pi0.V());
 
    // pi0 decay to 2gamma, in pi0 rest frame
    double gam1_p = CMmomentum (PI0_MASS, 0.0, 0.0);
    gamma1.polar(gam1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    gamma1.t(gam1_p);
    gamma2.set(gam1_p,threeVec (0.0, 0.0, 0.0) - gamma1.V());
    gamma2.t(gam1_p); 

    /*
     *  transform all 4-vectors back to lab frame
     */
    fourVec tmp;

    // boost gammas to pi0 rest frame
    tmp.set(pi0.t(),zeroVec - pi0.V());
    Boost.set(tmp);
    gamma1 = Boost * gamma1;
    gamma2 = Boost * gamma2;


    // boost from CM to target rest frame (lab)
    Boost.set (CMtarget);
    resonance *= Boost;
    pi0 *= Boost;
    proton *= Boost;
    gamma1 = Boost * gamma1;
    gamma2 = Boost * gamma2;

    /*
     *  generate vertices
     */
    threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
    //	threeVec decay;
    if (Print) {
      // std::cerr << "\n\n*** New Event\n";
      // std::cerr << "Beam:\n  " << beam;
      // std::cerr << "Beam in CM:\n  " << CMbeam;
      // std::cerr << "Target in CM:\n  " << CMtarget;
      // std::cerr << "Resonance\n";
      // std::cerr << "  Resonance mass: " << resonance_mass << "\n";
      // std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
      // std::cerr << "  Resonance:\n     " << resonance;
      // std::cerr << "  pi0\n";
      // std::cerr << "  pi0 mass: " << PI0_MASS << "\n";
      // std::cerr << "  pi0 momentum: " << pi0_p << "\n";
      // std::cerr << "  pi0:\n    " << pi0;
      // std::cerr << "  proton :\n  " << proton;
      // std::cerr << "vertices:\n";
      // std::cerr << "  prod: " << production;
    }

    if (debug) {

      std::cout << (proton + pi0).lenSq() << " ";
      std::cout << std::endl;

    }

    else {

      /*
       *  write event
       */	

      if (txt2part_style) {
	if (printGamma)
	  std::cout << "3" << std::endl;
	else
	  std::cout << "2" << std::endl;
	std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;
	pParticle_txt2part(Proton, production, proton);
	if (printGamma) {
	  pParticle_txt2part(Gamma, production, gamma1);
	  pParticle_txt2part(Gamma, production, gamma2);
	}
	else
	  pParticle_txt2part(Pi0, production, pi0);
      }
      else {
	if (printBeam)
	  pParticle(Gamma,production,beam);


	if (printBaryon)
	  pParticle(Proton,production,proton);
	if (printGamma) {
	  pParticle(Gamma,production,gamma1);
	  pParticle(Gamma,production,gamma2);
	}
	else {
	  pParticle(Pi0,production,pi0);
	}
      }
    }
    maxevents--;
  }
}







// resonance -> pi+ pi-
// isobar1 -> piplus + piminus2

void ppi0 (int argc, char *argv[])
{
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double masslow = 2 * PI_MASS,
    masshigh = 3.0,
    r2 = 0.0,
    t_min = 0.0,
    t_max,
    slope = 3.0,
    LorentzFactor = 0,
    expt_min = exp (-slope * t_min),
    expt_max,
    lfmax = 0,
    resonance_mass;
  fourVec beam,
    target,
    electron,
    resonance,
    proton,
    pi0;

  vector3_t 
    vbeam,
    pbeam;

  lorentzTransform Boost;

    for (int iarg = 1; iarg < argc; ++iarg) {
      char *ptr = argv[iarg];
      if (*ptr == '-') {
	ptr++;
	switch (*ptr) {
	case 'm':
	  ptr++;
	  maxevents = atoi (ptr);
	  std::cerr << "maxevents: " << maxevents << "\n";
	  break;
	case 'l':
	  ptr++;
	  lfevents = atoi (ptr);
	  std::cerr << "lfevents: " << lfevents << "\n";
	  break;
	case 'L':
	  ptr++;
	  masslow = atof (ptr);
	  break;
	case 'U':
	  ptr++;
	  masshigh = atof (ptr);
	  break;
	case 'p':
	  Print = 1;
	  break;
	case 'M':
	  break;

	case 'h':
	  UsageM2 (argv[0]);
	  return;
	default:
	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
	  break;
	}
      }
    }
  

  /*
   *-- beam and target in lab frame
   */


  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec (sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) +
                                                       pow((double) pbeam.z,2.0) + pow (BEAM_MASS, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));
    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();
   double n_pi_threshold = PROTON_MASS + PI_MASS;

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PI_MASS;

    if (masshigh < n_pi_threshold) {
      std::cerr << "p pi high mass below p pi threshold" << std::endl;
      exit(1);
    }
   if (masshigh > CMenergy - PI_MASS)
      masshigh = CMenergy - PI_MASS;
    if (masslow > CMenergy - PI_MASS) {
      std::cerr << "Not enough beam energy... exiting" << std::endl;
      exit(1);
    }

    t_max = pow (CMmomentum (2 * CMenergy, masslow, PROTON_MASS), 2.0);
    expt_max = exp (-slope * t_max);

    /*
     *-- generate the resonance
     */
    do {
      resonance_mass = randm (masslow, masshigh);
    } while (resonance_mass < n_pi_threshold);

    r2 = randm (0.0, 1.0);

    double resonance_p = CMmomentum (CMenergy, resonance_mass, PI_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double expt;
    double t;
    double tmin; 
    tmin = pow(resonance_mass,2.0) + pow(PROTON_MASS,2) - 2.0 * CMtarget.t() * resonance_E + 2.0 * ~CMtarget.V() * resonance_p;
    expt_min = exp (-slope * tmin);

    do {
      expt = randm (expt_min, expt_max);
      t = -log (expt) / slope;
      /* t -s really -t */
      costheta =- (CMtarget.t() * resonance_E
		  - 0.5 * (t + pow (PROTON_MASS, 2.0) + pow (resonance_mass, 2.0)))
	/ (~CMtarget.V() * resonance_p);
    } while (fabs (costheta) > 1.0);

    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t( resonance_E);



    /*
     *  now do decay in resonance rest frame
     */
    double pi0_p = CMmomentum (resonance_mass, PROTON_MASS, PI_MASS);

    pi0.polar(pi0_p,acos (randm (-0.9999, 0.9999)),randm(-M_PI,M_PI));
    pi0.t(sqrt (pi0.V().lenSq () + pow (PI_MASS, 2.0)));
    proton.V() = threeVec (0.0, 0.0, 0.0) - pi0.V();
    proton.t(sqrt (proton.V().lenSq () + pow (PROTON_MASS, 2.0)));

    /*
     *  compute lorentz factor
     */
    LorentzFactor = resonance_p * (pi0_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	//	fourVec tmpvec;
	/*
	 *  transform all 4-vectors back to lab frame
	 */
	fourVec tmp;
	threeVec zeroVec = threeVec (0.0, 0.0, 0.0);

	// boost from n_pi rest frame to CM
	tmp.V() = zeroVec - resonance.V();
	tmp.t(resonance.t());
	Boost.set (tmp);
	pi0 *= Boost;
	proton *= Boost;

	// boost from CM to target rest frame (lab)
	Boost.set (CMtarget);
	resonance *= Boost;
	pi0 *= Boost;
	proton *= Boost;

	/*
	 *  generate vertices
	 */
	// e852 values
	threeVec production = threeVec (0.0, 0.0, 0.0);
	//	threeVec decay;
	  electron = beam + target -proton - pi0;
	if (Print) {
	  // std::cerr << "\n\n*** New Event\n";
	  // std::cerr << "Beam:\n  " << beam;
	  // std::cerr << "Beam in CM:\n  " << CMbeam;
	  // std::cerr << "Target in CM:\n  " << CMtarget;
	  // std::cerr << "Resonance\n";
	  // std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  // std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  // std::cerr << "  expt: " << expt << "\n";
	  // std::cerr << "  t: " << t << "\n";
	  // std::cerr << " Qsq: "  << -(beam - electron).lenSq() << std::endl; 
	  // std::cerr << " t: " << (target - proton - pi0).lenSq() << std::endl;
	  // std::cerr << "  Resonance:\n     " << resonance;
	  // std::cerr << "pi0\n";
	  // std::cerr << "  pi0 mass: " << PI_MASS << "\n";
	  // std::cerr << "  pi0 momentum: " << pi0_p << "\n";
	  // std::cerr << "  pi0:\n    " << pi0;
	  // std::cerr << "  proton :\n  " << proton;
	  // std::cerr << "  electron :\n  " << electron;
	  // std::cerr << "vertices:\n";
	  // std::cerr << "  prod: " << production;
	  // std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cerr << pow (~(pi0 + proton), 2.0) << std::endl;
	}
	/*
	 *  write event
	 */
	  pParticle(Electron,production,electron);
	  pParticle(Proton,production,proton);
	maxevents--;
      }
    }
  }
}



/* ---------------- K K pi -------------------------------------*/

/*  K+ Ks Pi- */
void kpkspim(int argc, char *argv[])
{
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double  masslow = KCHARGED_MASS + KZERO_MASS + PI_MASS,
    masshigh = 0,
    t_max,   t_min = 0.0,
    slope = 3.5,  
    expt_min = exp (-slope * t_min),
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec
    beam,
    target,
    resonance,
    recoil,
    kplus,
    kshort,
    piminus,
    pi1, pi2,
    isobar1;
  lorentzTransform Boost;
  vector3_t vbeam,pbeam;
  float beamMass = BEAM_MASS;


  if (argc <= 1) {
    UsageM5(argv[0]);
    return;
  }
  else {
    for (int iarg = 1; iarg < argc; ++iarg) {
      char *ptr = argv[iarg];
      if (*ptr == '-') {
        ptr++;
        switch (*ptr) {
        case 'm':
          ptr++;
          maxevents = atoi (ptr);
          std::cerr << "maxevents: " << maxevents << "\n";
          break;
        case 'l':
          ptr++;
          lfevents = atoi (ptr);
          std::cerr << "lfevents: " << lfevents << "\n";
          break;
        case 'L':
          ptr++;
          masslow = atof (ptr);
          break;
        case 'U':
          ptr++;
          masshigh = atof (ptr);
          break;
         case 'p':
          Print = 1;
          break;
        case 'M':
          break;
	case 'b':
	  beamMass = atof(++ptr);
	  break;
        case 'h':
          UsageM5 (argv[0]);
          return;
        default:
	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
          break;
        }
      }
    }
  }


  double k_k_pi_threshold = KCHARGED_MASS + KZERO_MASS + PI_MASS;
  double k_pi_threshold =  KZERO_MASS + PI_MASS;
  if(masslow < k_k_pi_threshold)
    masslow = k_k_pi_threshold;


  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);


    /*
     *-- beam and target in lab frame
     */
/* beam mass changed to K- mass (kdb 11-18-96) */
    beam = fourVec (sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass,2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));

    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));
    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PROTON_MASS;
    if (masshigh < k_k_pi_threshold) {
      std::cerr << "k_k_pi high mass below k_k_pi_threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < k_k_pi_threshold);

    isobar1_mass = randm(k_pi_threshold, (resonance_mass - KCHARGED_MASS));

    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double expt;
    double t;  

    t_max = pow (CMmomentum (2 * CMenergy, masslow, PROTON_MASS), 2.0);   
    expt_max = exp (-slope * t_max);

    do {
      expt = randm (expt_min, expt_max);
      t = -log (expt) / slope;
      costheta = (CMbeam.t() * resonance_E
		  - 0.5 * (t + pow (PI_MASS, 2.0) + pow (resonance_mass, 2.0)))
	/ (~CMbeam.V() * resonance_p);
    } while (fabs (costheta) > 1.0);

    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);



    /*
     *-- recoil particle
     */
    recoil.V() = threeVec (0.0, 0.0, 0.0) - (resonance.V());
    recoil.t(sqrt (recoil.V().lenSq () + pow (PROTON_MASS, 2.0)));

    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, KCHARGED_MASS); 
    isobar1.polar(isobar1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI));
    isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    kplus.V() = threeVec (0.0, 0.0, 0.0) - isobar1.V();
    kplus.t(sqrt (kplus.V().lenSq () + pow (KCHARGED_MASS, 2.0)));

    /*
     *  now do decay in isobar1 rest frame
     */
  // k pi case
    double kshort_p = CMmomentum (isobar1_mass, KZERO_MASS,PI_MASS);
    kshort.polar(kshort_p,acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI));
    kshort.t(sqrt (kshort.V().lenSq () + pow (KCHARGED_MASS, 2.0)));
    piminus.V() = threeVec (0.0, 0.0, 0.0) - kshort.V();
    piminus.t(sqrt (piminus.V().lenSq () + pow (PI0_MASS, 2.0)));

         /* Kshort decay */


        double pi1_p = CMmomentum (KZERO_MASS, PI_MASS,PI_MASS);

	pi1.polar( pi1_p,acos (randm (-1.0, 1.0)), randm (-M_PI, M_PI));
        pi1.t(sqrt (pi1.V().lenSq () + pow (PI_MASS, 2.0)));
        pi2.V() = threeVec (0.0, 0.0, 0.0) - pi1.V();
        pi2.t(pi1.t());


    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar1_p * kshort_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
        /*
         *  transform all 4-vectors back to lab frame
         */
        fourVec tmp;
        threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
         /* Boost the pions  to Isobar  rest frame */
        tmp.V() = zeroVec - kshort.V();
        tmp.t(kshort.t());
        Boost.set (tmp);
        pi1 = Boost * pi1;
        pi2 = Boost * pi2;
 
        // boost to resonance rest frame

        tmp.V() = zeroVec - isobar1.V();
        tmp.t(isobar1.t());
        Boost.set (tmp);
        kshort = Boost * kshort;
        piminus = Boost * piminus;
        pi1 = Boost * pi1;
        pi2 = Boost * pi2;


        // boost from resonance rest frame to CM
        tmp.V() = zeroVec - resonance.V();
        tmp.t(resonance.t());
        Boost.set (tmp);
        isobar1 = Boost * isobar1;
        kshort = Boost * kshort;
        kplus = Boost * kplus;
        piminus = Boost * piminus;
        pi1 = Boost * pi1;
        pi2 = Boost * pi2;


        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        isobar1 = Boost * isobar1;
        kshort = Boost * kshort;
        kplus = Boost * kplus;
        piminus = Boost * piminus;
        pi1 = Boost * pi1;
        pi2 = Boost * pi2;

        /*
         *  generate vertices
         */
        // e852 values
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
    //    threeVec decay;



        if (Print) {
          // std::cerr << "\n\n*** New Event\n";
          // std::cerr << "Beam:\n  " << beam;
          // std::cerr << "Beam in CM:\n  " << CMbeam;
          // std::cerr << "Target in CM:\n  " << CMtarget;
          // std::cerr << "Resonance\n";
          // std::cerr << "  Resonance mass: " << resonance_mass << "\n";
          // std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  // std::cerr << "  t: " << t << "\n";
          // std::cerr << "  Resonance:\n     " << resonance;
          // std::cerr << "recoil: \n  " << recoil;
          // std::cerr << "isobar1\n";
          // std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
          // std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
          // std::cerr << "  isobar1:\n    " << isobar1;
          // std::cerr << "kshort :\n  " << kshort;
          // std::cerr << "pi-:\n  " << piminus;
          // std::cerr << "k+:\n  " << kplus;
          // std::cerr << "pi1:\n  " << pi1;
          // std::cerr << "pi2:\n  " << pi2;
          // std::cerr << "pi pi mass: " << ~(pi1 + pi2) << std::endl;
          // std::cerr << "vertices:\n";
          // std::cerr << "  prod: " << production;
          // std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
        }

        /*
         *  write event
         */

	pParticle(PiMinus,production,pi1);
	pParticle(PiPlus,production,pi2);
	pParticle(KPlus,production,kplus);
	pParticle(PiMinus,production,piminus);
	  
        maxevents--;
      }
    }
  }
}

/* ---------------- gamma p -> pi pi pi -------------------------------------*/

/*  pi+ pi+ pi- n */

void pipipi_gamma (int argc, char *argv[])
{
  int icount = 0;
  int Print = 0,
  maxevents = 9999999,
  nevents = 0,
  lfevents = 5000;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    slope = 10.0,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    piplus1,
    piminus,
    piplus2,
    isobar1;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = GAMMA_MASS;
  int printBaryon = 0;

  float tMin;

    for (int iarg = 1; iarg < argc; ++iarg) {
      char *ptr = argv[iarg];
      if (*ptr == '-') {
	ptr++;
	switch (*ptr) {
	case 'm':
	  ptr++;
	  maxevents = atoi (ptr);
	  std::cerr << "maxevents: " << maxevents << "\n";
	  break;
	case 'l':
	  ptr++;
	  lfevents = atoi (ptr);
	  std::cerr << "lfevents: " << lfevents << "\n";
	  break;
	case 'L':
	  ptr++;
	  masslow = atof(ptr);
	  break;
	case 'U':
	  ptr++;
	  masshigh = atof(ptr);
	  break;
	case 'p':
	  Print = 1;
	  break;
	case 'M':
	  break;
	case 'b':
	  beamMass = atof(++ptr);
	  break;
	case 'B':
	  printBaryon = 1;
	  break;
	case 'h':
	  UsageM7(argv[0]);
	  return;
	default:
	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
	  break;
	}
      }
    }


  double pi_pi_pi_threshold = 3 * PI_MASS;
  double pi_pi_threshold = 2 * PI_MASS;
  if(masslow < pi_pi_pi_threshold)
    masslow = pi_pi_pi_threshold;

 
   while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */

    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PROTON_MASS;
    if (masshigh < pi_pi_pi_threshold) {
      std::cerr << "pi_pi_pi high mass below pi_pi_pi_threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < pi_pi_pi_threshold);
    isobar1_mass = randm(pi_pi_threshold, (resonance_mass - PI_MASS));

    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,NEUTRON_MASS);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;
    icount = 0;
      if (TMax > 0.0) {
	do {
	  float r;
	r = randm(0.0,1.0);
	t = -log(r)/slope;
	}
	while((t < TLow) ||  (t > TMax));
      }
      else {


    icount = 0;
      do {
	if((icount++) > 100000) {
	  std::cerr << "exiting... too many counts" << std::endl;
	  exit(0);
	}
      

      t = randm(0., t_max);
      //     expt_max = t_max;
    } while( randm(0.,expt_max) > t*exp(-slope*t) );
      }
    costheta = 1. - 2.*t/t_max;
    //   costheta = 1.0;

    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, PI_MASS);
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

     isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    piplus1.set(sqrt(isobar1.V().lenSq() + pow (PI_MASS, 2.0)),zeroVec -isobar1.V());




    /*
     *  now do decay in isobar1 rest frame
     */
  // pi pi case
    double piminus_p = CMmomentum (isobar1_mass, PI_MASS,PI_MASS);

    piminus.polar( piminus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   piminus.polar( piminus_p,acos (1.0),0.0);
    piminus.t(sqrt (piminus.V().lenSq () + pow (PI_MASS, 2.0)));
    piplus2.set(sqrt (piminus.V().lenSq () + pow (PI_MASS, 2.0)),zeroVec - piminus.V());


    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar1_p * piminus_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
         /* transform all 4-vectors back to lab frame */
        fourVec tmp;

        // boost from pi_pi(rho) rest frame to pi_pi_pi(a2/M) rest frame

        tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        piminus = Boost * piminus;
        piplus2 = Boost * piplus2;

        // boost from pi_pi_pi rest frame to CM
        tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);

        isobar1 = Boost * isobar1;
        piminus = Boost * piminus;
        piplus2 = Boost * piplus2;
        piplus1 = Boost * piplus1;

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        isobar1 = Boost * isobar1;
        piminus = Boost * piminus;
        piplus2 = Boost * piplus2;
        piplus1 = Boost * piplus1;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

	if (Print) {
	   std::cerr << "\n\n*** New Event\n";
	   std::cerr << "Beam:\n  ";
	   beam.print();
	   std::cerr << "Beam in CM:\n  ";
	   CMbeam.print();
	   std::cerr << "Target in CM:\n  " ;
	   CMtarget.print();
	   std::cerr << "Resonance\n";
	   std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	   std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	   std::cerr << "  t: " << t << "\n";
	   std::cerr << "  Resonance:\n "    ;
	   resonance.print();
	   std::cerr << "recoil: \n  ";
	   recoil.print();
	   std::cerr << "isobar1\n";
	   std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	   std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	   std::cerr << "  isobar1:\n    ";
	   isobar1.print();
	   std::cerr << "pi+ :\n  ";
	   piminus.print();
	   std::cerr << "pi-2:\n  ";
	   piplus2.print();
	   std::cerr << "pi-1:\n  ";
	   piplus1.print();
	   //	   std::cerr << "vertices:\n";
	   //	   std::cerr << "  prod: " << production;
	   std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	   std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	  if (dalitz) {
	    std::cout << resonance_mass << " " << isobar1_mass*isobar1_mass  << " ";
	    std::cout << ~(piplus2 + piplus1 + piminus) << " ";
	    std::cout << pow (~(piplus2 + piminus), 2.0) << " "
	      << pow (~(piplus2 + piplus1), 2.0) << " "
		<< pow (~(piminus + piplus1), 2.0) << " "
		 << (piminus + recoil).lenSq() << " " 
		 << (piplus1 + recoil).lenSq() << " " 
		 << (piplus2 + recoil).lenSq() << " " ;
	    std::cout << masslow << " " << masshigh << " ";
	    std::cout << std::endl;
	  }

	  else if (txt2part_style){
	    std::cout << "4" << std::endl;
	  
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Neutron, production, recoil);
	  pParticle_txt2part(PiMinus, production, piminus);
	  pParticle_txt2part(PiPlus, production, piplus1);
	  pParticle_txt2part(PiPlus, production, piplus2);
	}

	  else if (gamp) {
	    std::cout << (printBaryon ? "5" : "4") << std::endl;
	    pParticleGamp(Gamma,beam);
	  if (printBaryon)
	    pParticleGamp(Neutron,recoil);
	  pParticleGamp(PiPlus,piplus1);
	  pParticleGamp(PiMinus,piminus);
	  pParticleGamp(PiPlus,piplus2);

	  }
	  else {
	/*
	 *  write event
	 */


	  if (printBaryon)
	    pParticle(Neutron,production,recoil);
	  pParticle(PiPlus,production,piplus1);
	  pParticle(PiMinus,production,piminus);
	  pParticle(PiPlus,production,piplus2);
	  }

	  maxevents--;
      }
    }
  }
}

/* ---------------- gamma p -> pi pi p -------------------------------------*/

/*  (pi+ pi-) p  */

void ppipi_gamma(int argc, char *argv[])
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double 
    masslow = 0.0,
    masshigh = 0., massHigh,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    piminus,
    piplus;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = GAMMA_MASS;
  int printBaryon = 0;
  int printGamma = 0;
  int init = 1;

  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break;
      case 'h':
	UsageM11(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 10.0;
  double pi_pi_threshold = 2 * PI_MASS;
  if(masslow < pi_pi_threshold)
    masslow = pi_pi_threshold;

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      massHigh = CMenergy - PROTON_MASS;
    else
      massHigh = masshigh;
    if (massHigh < pi_pi_threshold) {
      std::cerr << "pi_pi high mass below pi_pi_threshold" << std::endl;
      //     exit(1);
    }

    if (init) {
       std::cerr << masslow << " < m(pi+ pi-) < " << massHigh << std::endl;
       init = 0;
    }
       
 

    do {
      resonance_mass = randm(masslow, massHigh);
    } while (resonance_mass <pi_pi_threshold);
 
    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,NEUTRON_MASS);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {


      icount = 0;
      do {
	if((icount++) > 100000) {
	  std::cerr << "exiting... too many counts" << std::endl;
	  exit(0);
	}

	t = randm(0., t_max);
      } while( randm(0.,expt_max) > t*exp(-slope*t) );
      costheta = 1. - 2.*t/t_max;
    }


    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    // pi pi case
    double piplus_p = CMmomentum (resonance_mass, PI_MASS,PI_MASS);

    piplus.polar( piplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 

    piplus.t(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)));
    piminus.set(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)),zeroVec - piplus.V());
    
    /*
     *  compute lorentz factor
     */
    LorentzFactor = Isotropic ? piplus_p * resonance_p : piplus_p;
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

        // boost from pi_pi rest frame to CM
        tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);
	
	
        piplus = Boost * piplus;
        piminus = Boost * piminus;

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        piplus = Boost * piplus;
        piminus = Boost * piminus;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "pi+ :\n  ";
	  piplus.print();
	  std::cerr << "pi-:\n  ";
	  piminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (txt2part_style){
	  if (printGamma) {
	    std::cout << "5" << std::endl;
	  } 
	  else {
	    std::cout << "4" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Proton, production, recoil);
	  pParticle_txt2part(PiMinus, production, piminus);
	  pParticle_txt2part(PiPlus, production, piplus);
	}
	else if (gamp) {
	  std::cout << (printBaryon ? "4" : "3") << std::endl;
	  pParticleGamp(Gamma,beam);
	  if (printBaryon)
	    pParticleGamp(Proton,recoil);
	  pParticleGamp(PiMinus,piminus);
	  pParticleGamp(PiPlus,piplus);
	}
	else if (dalitz) {
	  std::cout << (recoil + piminus).lenSq() << " "  << (recoil + piplus).lenSq() << " " << (piplus + piminus).lenSq() << std::endl;
	}
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Gamma,production,beam);


	  if (printBaryon)
	    pParticle(Proton,production,recoil);
	  pParticle(PiPlus,production,piminus);
	  pParticle(PiMinus,production,piplus);
	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

/*----------------End of 3 pi-----------------------------------*/

/* ---------------- gamma p -> e+ e- p -------------------------------------*/

/*  (e+ e-) p  */

void pepem_gamma(int argc, char *argv[])
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double 
    masslow = 2. * ELEC_MASS,
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    eminus,
    eplus;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = GAMMA_MASS;
  int printBaryon = 0;
  int printGamma = 0;

  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break;
      case 'h':
	UsageM15(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 10.0;
  double e_e_threshold = 2. * ELEC_MASS;
  if(masslow < e_e_threshold)
    masslow = e_e_threshold;
  if (!Isotropic)
    std::cerr << "Generate events with t-slope " << slope << std::endl;


  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pbeam.x * pbeam.x + pbeam.y * pbeam.y + pbeam.z *  pbeam.z + beamMass * beamMass),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PROTON_MASS;
    if (masshigh < e_e_threshold) {
      std::cerr << "e+_e- high mass below e_e_threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass <e_e_threshold);
 
    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (resonance_mass * resonance_mass + resonance_p * resonance_p);
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,NEUTRON_MASS);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {


    icount = 0;
    do {
      if((icount++) > 100000) {
	std::cerr << "exiting... too many counts" << std::endl;
	exit(0);
      }

      t = randm(0., t_max);
    } while( randm(0.,expt_max) > exp(-slope*t) );
    costheta = 1. - 2.*t/t_max;
    }


    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  PROTON_MASS * PROTON_MASS),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    // e+ e- case
    double eplus_p = CMmomentum (resonance_mass, ELEC_MASS,ELEC_MASS);

    eplus.polar( eplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 

    eplus.t(sqrt (eplus.V().lenSq () + ELEC_MASS * ELEC_MASS ));
    eminus.set(sqrt (eplus.V().lenSq () + ELEC_MASS * ELEC_MASS),zeroVec - eplus.V());
    
    /*
     *  compute lorentz factor
     */
    LorentzFactor = eplus_p;
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

        // boost from pi_pi rest frame to CM
        tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);
	
	
        eplus = Boost * eplus;
        eminus = Boost * eminus;

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        eplus = Boost * eplus;
        eminus = Boost * eminus;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "e+ :\n  ";
	  eplus.print();
	  std::cerr << "e-:\n  ";
	  eminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << resonance_mass << " " << isobar1_mass*isobar1_mass  << " ";
	  std::cout << ~(eminus + eplus) << " ";
	    std::cout   << pow (~(eplus + eminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << (eminus + eplus + recoil).lenSq() << " ";
	  std::cout << std::endl;
	} 
	else if (txt2part_style){
	    std::cout << "3" << std::endl;
	  
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Proton, production, recoil);
	  pParticle_txt2part(Electron, production, eminus);
	  pParticle_txt2part(Positron, production, eplus);
	}

	else if (gamp) {
	    std::cout << (printBaryon ? "4" : "3") << std::endl;
	    pParticleGamp(Gamma,beam);
	  if (printBaryon)
	    pParticleGamp(Proton,recoil);
	  pParticleGamp(Electron,eminus);
	  pParticleGamp(Positron,eplus);
	}
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Gamma,production,beam);


	  if (printBaryon)
	    pParticle(Proton,production,recoil);
	  pParticle(Electron,production,eminus);
	  pParticle(Positron,production,eplus);
	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

/*---------------- End of p e+ e- -----------------------------------*/

/* ---------------- gamma p -> pi pi p -------------------------------------*/

/*  (pi+ pi-) p  */

void ppipiX_gamma(int argc, char *argv[])
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    piminus,
    piplus;
  lorentzTransform Boost;

  char mode = '-';
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = GAMMA_MASS;
  int printBaryon = 0;

  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'c':
	mode = *(++ptr);
	break;
      case 'h':
	UsageM11(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 10.0;
  double p_pi_threshold = PI_MASS + PROTON_MASS;
  if(masslow < p_pi_threshold)
    masslow = p_pi_threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - PI_MASS;
    if (masshigh < p_pi_threshold) {
      std::cerr << "p_pi high mass below p_pi_threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass <p_pi_threshold);
 
    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, PI_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,NEUTRON_MASS);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {


    icount = 0;
    do {
      if((icount++) > 100000) {
	std::cerr << "exiting... too many counts" << std::endl;
	exit(0);
      }

      t = randm(0., t_max);
    } while( randm(0.,expt_max) > t*exp(-slope*t) );
    costheta = -(1. - 2.*t/t_max);
    }


    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle (the pi-)
     */
    piminus.set(sqrt(resonance.V().lenSq() +  pow (PI_MASS, 2.0)),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    // pi pi case
    double piplus_p = CMmomentum (resonance_mass, PROTON_MASS,PI_MASS);

    piplus.polar( piplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 

    piplus.t(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)));
    recoil.set(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)),zeroVec - piplus.V());
    
    /*
     *  compute lorentz factor
     */
    LorentzFactor = piplus_p;
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

        // boost from p_pi rest frame to CM
        tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);
	
	
        piplus = Boost * piplus;
        recoil = Boost * recoil;

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        piplus = Boost * piplus;
        piminus = Boost * piminus;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

	if (mode == '+') {
	  tmp = piminus;
	  piminus = piplus;
	  piplus = tmp;
	}
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "pi+ :\n  ";
	  piplus.print();
	  std::cerr << "pi-:\n  ";
	  piminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << resonance_mass << " " << isobar1_mass*isobar1_mass  << " ";
	  std::cout << ~(piminus + piplus) << " ";
	    std::cout   << pow (~(piplus + piminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << (piminus + piplus + recoil).lenSq() << " ";
	  std::cout << std::endl;
	} 
	else if (txt2part_style){
	  std::cout << "3" << std::endl;
	  std::cout << EbeamZ <<  " " << beam.t() << " " << production.z()/SPEED_OF_LIGHT << std::endl;
	  pParticle_txt2part(Proton, production, recoil);
	  pParticle_txt2part(PiMinus, production, piminus);
	  pParticle_txt2part(PiPlus, production, piplus);
	}
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Gamma,production,beam);


	  if (printBaryon)
	    pParticle(Proton,production,recoil);
	  pParticle(PiPlus,production,piminus);
	  pParticle(PiMinus,production,piplus);
	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

/*----------------End of 3 pi-----------------------------------*/

/* ---------------- gamma p -> omega pi pi  -------------------------------------*/

/*  omega pi pip */

void omegapipi (int argc, char *argv[],Particle_t Beam,Particle_t Pi1,Particle_t Pi2)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  int PrintGammas = 0;
  Particle_t Baryon;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    omeg,
    pi1,
    pi2,
    piplus,
    piminus,
    pizero,
    gamma1,
    gamma2;
  fourVec isobar1,isobar2;
  // isobar1 is pi pi isobar off of omega pi pi system
  // isobar2 is pi+ pi- isobar from omega decay
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = GAMMA_MASS;
  int printBaryon = 0;
  int printGamma = 0;

  double pi1mass  = Mass(Pi1),pi2mass = Mass(Pi2),baryonmass = Mass(Baryon);

  float tMin;

  int debug = 0;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'D':
	debug = 1;
	break;
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break;
      case 'h':
	UsageM16(argv[0]);
	return;
      default:
	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 4.0;
  double threshold = OMEGA_MASS + pi1mass + pi2mass;
  double pi_pi_threshold = pi1mass + pi2mass;

  int qtot = Q(Beam) + Q(Proton) - Q(Pi1) - Q(Pi2);

  switch (qtot) {
  case 0:
    Baryon = Neutron;
    break;
  case 1:
    Baryon = Proton;
    break;
  default:
    std::cerr << "illegal charge combination:\ttotal charge  = " << qtot << std::endl;
    exit(0);
  }


  beamMass = Mass(Beam);
  baryonmass = Mass(Baryon);

  // set up the recoil

  
  if(masslow < threshold)
    masslow = threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - baryonmass;
    if (masshigh < threshold) {
      std::cerr << "resonance  high mass below threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < threshold);
    isobar1_mass = randm(pi_pi_threshold, (resonance_mass - OMEGA_MASS));

    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, baryonmass);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {


    icount = 0;
    do {
      if((icount++) > 100000) {
	std::cerr << "exiting... too many counts" << std::endl;
	exit(0);
      }

      t = randm(0., t_max);
      //     expt_max = t_max;
    } while( randm(0.,expt_max) > t*exp(-slope*t) );
    costheta = 1. - 2.*t/t_max;
    //     costheta = 1.0;
    }


    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());



    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, OMEGA_MASS);
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

    isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    omeg.set(sqrt(isobar1.V().lenSq() + pow (OMEGA_MASS, 2.0)),zeroVec -isobar1.V());

    /*
     *  now do decay in isobar1 rest frame
     */
    // pi pi case
    double pi1_p = CMmomentum (isobar1_mass, pi1mass,pi2mass);

    pi1.polar( pi1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   piplus.polar( piplus_p,acos (1.0),0.0);
    pi1.t(sqrt (pi1.V().lenSq () + pow (pi1mass, 2.0)));
    pi2.set(sqrt (pi1.V().lenSq () + pow (pi2mass, 2.0)),zeroVec - pi1.V());


    // code to do omega decay
    double himass = OMEGA_MASS - PI0_MASS;
    double lomass = 2 * PI_MASS;
    double isobar2_mass = randm(lomass,himass);
    double isobar2_p = CMmomentum(OMEGA_MASS,isobar2_mass,PI0_MASS);
    isobar2.polar( isobar2_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI));
    isobar2.t(sqrt (isobar2.V().lenSq () + pow (isobar2_mass, 2.0))); 
    pizero.set(sqrt(isobar2.V().lenSq() + pow(PI0_MASS,2.0)),zeroVec - isobar2.V());

 
    // isobar2 decay
    double piplus_p = CMmomentum(isobar2_mass,PI_MASS,PI_MASS);

    piplus.polar( piplus_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    piplus.t(sqrt(piplus.V().lenSq() + pow(PI_MASS,2.0)));
    piminus.set(sqrt(piplus.V().lenSq() + pow(PI_MASS,2.0)),zeroVec - piplus.V());
	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi+ (isobar2) : "; pParticleGamp(PiPlus,piplus);
	  std::cout << "pi- (isobar2) : "; pParticleGamp(PiMinus,piminus);
	  std::cout << "mass " << ~(piminus + piplus) << " " << ~isobar2 << " " << isobar2_mass << " " << piplus_p << " " << piplus.r() << std::endl;
	  std::cout << std::endl << std::endl;
	}



    
    // pi0 decay to 2 gamma, in pi0 rest frame
    double gam1_p = CMmomentum (PI0_MASS, 0.0, 0.0);
    gamma1.polar(gam1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    gamma1.t(gam1_p);
    gamma2.set(gam1_p,threeVec (0.0, 0.0, 0.0) - gamma1.V());
    gamma2.t(gam1_p);
    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar1_p * isobar2_p * pi1_p * piplus_p * gam1_p );
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

	// boost gammas to omega  rest frame
	tmp.set(pizero.t(),zeroVec - pizero.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;

	// boost pi+ pi- to omega rest frame
	tmp.set(isobar2.t(),zeroVec - isobar2.V());
	Boost.set(tmp);



	piplus = Boost * piplus;
	piminus = Boost * piminus;

	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi+ (isobar2) : "; pParticleGamp(PiPlus,piplus);
	  std::cout << "pi- (isobar2) : "; pParticleGamp(PiMinus,piminus);
	  std::cout << "mass " << ~(piminus + piplus) << " " << ~isobar2 << " " << isobar2_mass << std::endl;
	  std::cout << std::endl << std::endl;
	}



	// boost omega progeny
	tmp.set(omeg.t(),zeroVec - omeg.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	piplus = Boost * piplus;
	piminus = Boost * piminus;
	pizero = Boost * pizero;
	isobar2 = Boost * isobar2;
	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi1 : "; pParticleGamp(Pi1,pi1);
	  std::cout << "pi2 : "; pParticleGamp(Pi2,pi2);
	  std::cout << "mass " << ~(pi1 + pi2) << std::endl;
	  std::cout << std::endl << std::endl;
	}

	// boost progeny of isobar1 
	tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        pi1 = Boost * pi1;
	pi2 = Boost * pi2;

	if (debug) {
	  std::cout << "\n";
	  std::cout << "isobar1: " << ~isobar1 << " "; pParticleGamp(Rho0,isobar1);
	  std::cout << "pi1 : "; pParticleGamp(Pi1,pi1);
	  std::cout << "pi2 : "; pParticleGamp(Pi2,pi2);
	  std::cout << std::endl << std::endl;
	}

        // boost  everyone to resonance   frame

        tmp.set(resonance.t(),zeroVec - resonance.V());
        Boost.set (tmp);
        piplus = Boost * piplus;
	piminus = Boost * piminus;
        pizero = Boost * pizero;
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	pi1 = Boost * pi1;
	pi2 = Boost * pi2; 
	isobar2 = Boost * isobar2;
	omeg = Boost * omeg;
	isobar1 = Boost * isobar1;



        // boost from CM to target rest frame (lab)
        tmp.set (CMtarget.t(),zeroVec - CMtarget.V());
	Boost.set(CMtarget);
	CMtarget = Boost * CMtarget;
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        piplus = Boost * piplus;
	piminus = Boost * piminus;
        pizero = Boost * pizero;
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	pi1 = Boost * pi1;
	pi2 = Boost * pi2;
	omeg = Boost * omeg;
	isobar1 = Boost * isobar1;
	isobar2 = Boost * isobar2;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "isobar1\n";
	  std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	  std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	  std::cerr << "  isobar1:\n    ";
	  isobar1.print();
	  std::cerr << "pi+ :\n  ";
	  piplus.print();
	  std::cerr << "pi-2:\n  ";
	  pizero.print();
	  std::cerr << "pi-1:\n  ";
	  piminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << resonance_mass << " " << isobar1_mass*isobar1_mass  << " ";
	  std::cout << ~(pizero + piminus + piplus) << " ";
	  std::cout << pow (~(pizero + piplus), 2.0) << " "
	       << pow (~(pizero + piminus), 2.0) << " "
	       << pow (~(piplus + piminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << (pizero + piminus + piplus + recoil).lenSq() << " ";
	  std::cout << std::endl;
	} 
	else if (txt2part_style){
	  if (printGamma) {
	    std::cout << "7" << std::endl;
	  } 
	  else {
	    std::cout << "6" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Baryon, production, recoil);
	  pParticle_txt2part(Pi1, production, pi1);
	  pParticle_txt2part(Pi2, production, pi2);
	  pParticle_txt2part(PiMinus,production,piminus);
	  pParticle_txt2part(PiPlus,production,piplus);
	  if (printGamma){
	    pParticle_txt2part(Gamma, production, gamma1);
	    pParticle_txt2part(Gamma, production, gamma2);
	  }
	  else{
	    pParticle_txt2part(Pi0,production,pizero);
	  }
	}
	else if (gamp){
	  if (printGamma) {
	    std::cout << "8" << std::endl;
	  } 
	  else {
	    std::cout << "7" << std::endl;
	  }
	  pParticleGamp(Beam,beam);
	  pParticleGamp(Baryon,recoil);
	  pParticleGamp(Pi1, pi1);
	  pParticleGamp(Pi2, pi2);
	  pParticleGamp(PiMinus,piminus);
	  pParticleGamp(PiPlus,piplus);
	  if (printGamma){
	    pParticleGamp(Gamma, gamma1);
	    pParticleGamp(Gamma, gamma2);
	  }
	  else{
	    pParticleGamp(Pi0,pizero);
	  }
	}
	else if (debug) {
	  std::cout << "5" << std::endl;
	  std::cout << "beam "; pParticleGamp(Beam,beam);
	  std::cout << "recoil "; pParticleGamp(Baryon,recoil);
	  std::cout << "Resonance "; pParticleGamp(Rho0,resonance);
	  std::cout << "omega "; pParticleGamp(omega,omeg);
	  std::cout << "isobar1 "; pParticleGamp(Rho0,isobar1);
	  std::cout << "isobar2 "; pParticleGamp(Rho0,isobar2);
	  std::cout << "pi1 (isobar1) "; pParticleGamp(Pi1,pi1);
	  std::cout << "pi2 (isobar1) "; pParticleGamp(Pi2,pi2);
	  std::cout << "pizero (omega)"; pParticleGamp(Pi0,pizero);
	  std::cout << "piplus (omega) "; pParticleGamp(PiPlus,piplus);
	  std::cout << "piminus (omega) "; pParticleGamp(PiMinus,piminus);
	}
	  
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Beam,production,beam);


	  if (printBaryon)
	    pParticle(Baryon,production,recoil);
	  pParticle(Pi1,production,pi1);
	  pParticle(Pi2,production,pi2);
	  pParticle(PiPlus,production,piminus);
	  pParticle(PiMinus,production,piplus);
	  if (printGamma) {
	    pParticle(Gamma,production,gamma1);
	    pParticle(Gamma,production,gamma2);
	  }
	  else {
	    pParticle(Pi0,production,pizero);
	  }
	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

/* ---------------- gamma p -> omega pi pi  -------------------------------------*/

/*  omega pi pip */

void omegapipi0(int argc, char *argv[],Particle_t Beam,Particle_t Pi1,Particle_t Pi2)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  int PrintGammas = 0;
  Particle_t Baryon;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    omeg,
    pi1,
    pi2,
    piplus,
    piminus,
    pizero,
    gamma1,
    gamma2,
    gamma3,
    gamma4;
  fourVec isobar1,isobar2;
  // isobar1 is pi pi isobar off of omega pi pi system
  // isobar2 is pi+ pi- isobar from omega decay
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = GAMMA_MASS;
  int printBaryon = 0;
  int printGamma = 0;

  double pi1mass  = Mass(Pi1),pi2mass = Mass(Pi2),baryonmass = Mass(Baryon);

  float tMin;

  int debug = 0;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'D':
	debug = 1;
	break;
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break;
      case 'h':
	UsageM16(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 4.0;
  double threshold = OMEGA_MASS + pi1mass + pi2mass;
  double pi_pi_threshold = pi1mass + pi2mass;

  int qtot = Q(Beam) + Q(Proton) - Q(Pi1) - Q(Pi2);

  switch (qtot) {
  case 0:
    Baryon = Neutron;
    break;
  case 1:
    Baryon = Proton;
    break;
  default:
    std::cerr << "illegal charge combination:\ttotal charge  = " << qtot << std::endl;
    exit(0);
  }


  beamMass = Mass(Beam);
  baryonmass = Mass(Baryon);

  // set up the recoil

  
  if(masslow < threshold)
    masslow = threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - baryonmass;
    if (masshigh < threshold) {
      std::cerr << "resonance  high mass below threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < threshold);
    isobar1_mass = randm(pi_pi_threshold, (resonance_mass - OMEGA_MASS));

    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, baryonmass);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {


      icount = 0;
      do {
	if((icount++) > 100000) {
	  std::cerr << "exiting... too many counts" << std::endl;
	  exit(0);
	}

	t = randm(0., t_max);
	//     expt_max = t_max;
      } while( randm(0.,expt_max) > t*exp(-slope*t) );
      costheta = 1. - 2.*t/t_max;
      //     costheta = 1.0;
    }


    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());



    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, OMEGA_MASS);
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

    isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    omeg.set(sqrt(isobar1.V().lenSq() + pow (OMEGA_MASS, 2.0)),zeroVec -isobar1.V());

    /*
     *  now do decay in isobar1 rest frame
     */
    // pi pi case
    double pi1_p = CMmomentum (isobar1_mass, pi1mass,pi2mass);

    pi1.polar( pi1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   piplus.polar( piplus_p,acos (1.0),0.0);
    pi1.t(sqrt (pi1.V().lenSq () + pow (pi1mass, 2.0)));
    pi2.set(sqrt (pi1.V().lenSq () + pow (pi2mass, 2.0)),zeroVec - pi1.V());


    // code to do omega decay
    double himass = OMEGA_MASS - PI0_MASS;
    double lomass = 2 * PI_MASS;
    double isobar2_mass = randm(lomass,himass);
    double isobar2_p = CMmomentum(OMEGA_MASS,isobar2_mass,PI0_MASS);
    isobar2.polar( isobar2_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI));
    isobar2.t(sqrt (isobar2.V().lenSq () + pow (isobar2_mass, 2.0))); 
    pizero.set(sqrt(isobar2.V().lenSq() + pow(PI0_MASS,2.0)),zeroVec - isobar2.V());

 
    // isobar2 decay
    double piplus_p = CMmomentum(isobar2_mass,PI_MASS,PI_MASS);

    piplus.polar( piplus_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    piplus.t(sqrt(piplus.V().lenSq() + pow(PI_MASS,2.0)));
    piminus.set(sqrt(piplus.V().lenSq() + pow(PI_MASS,2.0)),zeroVec - piplus.V());
    if (debug) {
      std::cout << "\n";
      std::cout << "pi+ (isobar2) : "; pParticleGamp(PiPlus,piplus);
      std::cout << "pi- (isobar2) : "; pParticleGamp(PiMinus,piminus);
      std::cout << "mass " << ~(piminus + piplus) << " " << ~isobar2 << " " << isobar2_mass << " " << piplus_p << " " << piplus.r() << std::endl;
      std::cout << std::endl << std::endl;
    }



    
    // pi0 decay to 2 gamma, in pi0 rest frame pi0
    double gam1_p = CMmomentum (PI0_MASS, 0.0, 0.0);
    gamma1.polar(gam1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    gamma1.t(gam1_p);
    gamma2.set(gam1_p,threeVec (0.0, 0.0, 0.0) - gamma1.V());
    gamma2.t(gam1_p);

    gamma3.polar(gam1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    gamma3.t(gam1_p);
    gamma4.set(gam1_p,threeVec (0.0, 0.0, 0.0) - gamma3.V());
    gamma4.t(gam1_p);
    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar1_p * isobar2_p * pi1_p * piplus_p * gam1_p * gam1_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

	// boost gammas to from first pi0 decay
	tmp.set(pi2.t(),zeroVec - pi2.V());
	Boost.set(tmp);
	gamma3 = Boost * gamma3;
	gamma4 = Boost * gamma4;

	// boost gammas to omega  rest frame
	tmp.set(pizero.t(),zeroVec - pizero.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;

	// boost pi+ pi- to omega rest frame
	tmp.set(isobar2.t(),zeroVec - isobar2.V());
	Boost.set(tmp);



	piplus = Boost * piplus;
	piminus = Boost * piminus;

	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi+ (isobar2) : "; pParticleGamp(PiPlus,piplus);
	  std::cout << "pi- (isobar2) : "; pParticleGamp(PiMinus,piminus);
	  std::cout << "mass " << ~(piminus + piplus) << " " << ~isobar2 << " " << isobar2_mass << std::endl;
	  std::cout << std::endl << std::endl;
	}



	// boost omega progeny
	tmp.set(omeg.t(),zeroVec - omeg.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	piplus = Boost * piplus;
	piminus = Boost * piminus;
	pizero = Boost * pizero;
	isobar2 = Boost * isobar2;
	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi1 : "; pParticleGamp(Pi1,pi1);
	  std::cout << "pi2 : "; pParticleGamp(Pi2,pi2);
	  std::cout << "mass " << ~(pi1 + pi2) << std::endl;
	  std::cout << std::endl << std::endl;
	}

	// boost progeny of isobar1 
	tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        pi1 = Boost * pi1;
	pi2 = Boost * pi2;
	gamma3 = Boost * gamma3;
	gamma4 = Boost * gamma4;

	if (debug) {
	  std::cout << "\n";
	  std::cout << "isobar1: " << ~isobar1 << " "; pParticleGamp(Rho0,isobar1);
	  std::cout << "pi1 : "; pParticleGamp(Pi1,pi1);
	  std::cout << "pi2 : "; pParticleGamp(Pi2,pi2);
	  std::cout << std::endl << std::endl;
	}

        // boost  everyone to resonance   frame

        tmp.set(resonance.t(),zeroVec - resonance.V());
        Boost.set (tmp);
        piplus = Boost * piplus;
	piminus = Boost * piminus;
        pizero = Boost * pizero;
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	pi1 = Boost * pi1;
	pi2 = Boost * pi2; 
	isobar2 = Boost * isobar2;
	omeg = Boost * omeg;
	isobar1 = Boost * isobar1;
	gamma3 = Boost * gamma3;
	gamma4 = Boost * gamma4;



        // boost from CM to target rest frame (lab)
        tmp.set (CMtarget.t(),zeroVec - CMtarget.V());
	Boost.set(CMtarget);
	CMtarget = Boost * CMtarget;
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        piplus = Boost * piplus;
	piminus = Boost * piminus;
        pizero = Boost * pizero;
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	pi1 = Boost * pi1;
	pi2 = Boost * pi2;
	omeg = Boost * omeg;
	isobar1 = Boost * isobar1;
	isobar2 = Boost * isobar2;
	gamma3 = Boost * gamma3;
	gamma4 = Boost * gamma4;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "isobar1\n";
	  std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	  std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	  std::cerr << "  isobar1:\n    ";
	  isobar1.print();
	  std::cerr << "pi+ :\n  ";
	  piplus.print();
	  std::cerr << "pi-2:\n  ";
	  pizero.print();
	  std::cerr << "pi-1:\n  ";
	  piminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << resonance_mass << " " << isobar1_mass*isobar1_mass  << " ";
	  std::cout << ~(pizero + piminus + piplus) << " ";
	  std::cout << pow (~(pizero + piplus), 2.0) << " "
	       << pow (~(pizero + piminus), 2.0) << " "
	       << pow (~(piplus + piminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << (pizero + piminus + piplus + recoil).lenSq() << " ";
	  std::cout << std::endl;
	} 
	else if (txt2part_style){
	  if (printGamma) {
	    std::cout << "8" << std::endl;
	  } 
	  else {
	    std::cout << "6" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Baryon, production, recoil);
	  pParticle_txt2part(Pi1, production, pi1);
	  pParticle_txt2part(PiMinus,production,piminus);
	  pParticle_txt2part(PiPlus,production,piplus);
	  if (printGamma){
	    pParticle_txt2part(Gamma, production, gamma1);
	    pParticle_txt2part(Gamma, production, gamma2);
	    pParticle_txt2part(Gamma, production, gamma3);
	    pParticle_txt2part(Gamma, production, gamma4);
	  }
	  else{
	    pParticle_txt2part(Pi2, production, pi2);
	    pParticle_txt2part(Pi0,production,pizero);
	  }
	}
	else if (gamp){
	  if (printAll) {
	    std::cout << "11 " << std::endl;
	    pParticleGamp(Beam,beam);
	    pParticleGamp(Baryon,recoil);
	    pParticleGamp(Pi1, pi1);
	    pParticleGamp(Pi2, pi2);
	    pParticleGamp(Pi0,pizero);
	    pParticleGamp(PiMinus,piminus);
	    pParticleGamp(PiPlus,piplus);	
	    pParticleGamp(Gamma, gamma1);
	    pParticleGamp(Gamma, gamma2);
	    pParticleGamp(Gamma, gamma3);
	    pParticleGamp(Gamma, gamma4);
	  }
	  else {
	    
	    if (printGamma) {
	      std::cout << "9" << std::endl;
	    } 
	    else {
	      std::cout << "7" << std::endl;
	    }
	    pParticleGamp(Beam,beam);
	    pParticleGamp(Baryon,recoil);
	    pParticleGamp(Pi1, pi1);
	    pParticleGamp(PiMinus,piminus);
	    pParticleGamp(PiPlus,piplus);
	    if (printGamma){
	      pParticleGamp(Gamma, gamma1);
	      pParticleGamp(Gamma, gamma2);
	      pParticleGamp(Gamma, gamma3);
	      pParticleGamp(Gamma, gamma4);
	    }
	    else{
	      pParticleGamp(Pi2, pi2);
	      pParticleGamp(Pi0,pizero);
	    }
	  }
	}
	else if (debug) {
	  std::cout << "5" << std::endl;
	  std::cout << "beam "; pParticleGamp(Beam,beam);
	  std::cout << "recoil "; pParticleGamp(Baryon,recoil);
	  std::cout << "Resonance "; pParticleGamp(Rho0,resonance);
	  std::cout << "omega "; pParticleGamp(omega,omeg);
	  std::cout << "isobar1 "; pParticleGamp(Rho0,isobar1);
	  std::cout << "isobar2 "; pParticleGamp(Rho0,isobar2);
	  std::cout << "pi1 (isobar1) "; pParticleGamp(Pi1,pi1);
	  std::cout << "pi2 (isobar1) "; pParticleGamp(Pi2,pi2);
	  std::cout << "pizero (omega)"; pParticleGamp(Pi0,pizero);
	  std::cout << "piplus (omega) "; pParticleGamp(PiPlus,piplus);
	  std::cout << "piminus (omega) "; pParticleGamp(PiMinus,piminus);
	}
	  
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Beam,production,beam);


	  if (printBaryon)
	    pParticle(Baryon,production,recoil);
	  pParticle(Pi1,production,pi1);
	  pParticle(Pi2,production,pi2);
	  pParticle(PiPlus,production,piminus);
	  pParticle(PiMinus,production,piplus);
	  if (printGamma) {
	    pParticle(Gamma,production,gamma1);
	    pParticle(Gamma,production,gamma2);
	  }
	  else {
	    pParticle(Pi0,production,pizero);
	  }
	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}


/*  omega pi pip */

void omegaphi (int argc, char *argv[],Particle_t Beam)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  int PrintGammas = 0;
  Particle_t Baryon;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    omeg,
          phi,
    Kp,
    Km,
          piplus,
          piminus,
    pizero,
    gamma1,
    gamma2;
  fourVec isobar2;
  // isobar1 is pi pi isobar off of omega pi pi system
  // isobar2 is pi+ pi- isobar from omega decay
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = GAMMA_MASS;
  int printBaryon = 0;
  int printGamma = 0;


  float tMin;

  int debug = 0;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'D':
	debug = 1;
	break;
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break;
      case 'h':
	UsageM16(argv[0]);
	return;
      default:
	  std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 4.0;
  double threshold = OMEGA_MASS + PHI_MASS;
  

  int qtot = Q(Beam) + Q(Proton);

  switch (qtot) {
  case 0:
    Baryon = Neutron;
    break;
  case 1:
    Baryon = Proton;
    break;
  default:
    std::cerr << "illegal charge combination:\ttotal charge  = " << qtot << std::endl;
    exit(0);
  }


  beamMass = Mass(Beam);
  double baryonmass = Mass(Baryon);

  // set up the recoil

  
  if(masslow < threshold)
    masslow = threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - baryonmass;
    if (masshigh < threshold) {
      std::cerr << "resonance  high mass below threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < threshold);

    double beam_p      = CMmomentum (CMenergy, BEAM_MASS, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, baryonmass);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {

        double tp = tprimeExp(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass,slope);
        costheta = 1. - 2.*t/t_max;
        costheta = cosThetaExp(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass,slope);
    //     costheta = 1.0;
    }

    
    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());



    /*
     *  now do decay in resonance rest frame
     */
    double phi_p = CMmomentum (resonance_mass, PHI_MASS, OMEGA_MASS);
     
    phi.polar(phi_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

    phi.t(sqrt (phi.V().lenSq () + pow (PHI_MASS, 2.0)));
    
    omeg.set(sqrt(phi.V().lenSq() + pow (OMEGA_MASS, 2.0)),zeroVec - phi.V());

    /*
     *  now do decay in phi rest frame
     */
    // K+ K- case
    double Kp_p = CMmomentum (PHI_MASS, KCHARGED_MASS,KCHARGED_MASS);

    Kp.polar( Kp_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    
    Kp.t(sqrt (Kp.V().lenSq () + pow (KCHARGED_MASS, 2.0)));
    Km.set(sqrt (Kp.V().lenSq () + pow (KCHARGED_MASS, 2.0)),zeroVec - Kp.V());
    
   


    // code to do omega decay
    double himass = OMEGA_MASS - PI0_MASS;
    double lomass = 2 * PI_MASS;
    double isobar2_mass = randm(lomass,himass);
    double isobar2_p = CMmomentum(OMEGA_MASS,isobar2_mass,PI0_MASS);
    isobar2.polar( isobar2_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI));
    isobar2.t(sqrt (isobar2.V().lenSq () + pow (isobar2_mass, 2.0))); 
    pizero.set(sqrt(isobar2.V().lenSq() + pow(PI0_MASS,2.0)),zeroVec - isobar2.V());

 
    // isobar2 decay
    double piplus_p = CMmomentum(isobar2_mass,PI_MASS,PI_MASS);

    piplus.polar( piplus_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    piplus.t(sqrt(piplus.V().lenSq() + pow(PI_MASS,2.0)));
    piminus.set(sqrt(piplus.V().lenSq() + pow(PI_MASS,2.0)),zeroVec - piplus.V());
	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi+ (isobar2) : "; pParticleGamp(PiPlus,piplus);
	  std::cout << "pi- (isobar2) : "; pParticleGamp(PiMinus,piminus);
	  std::cout << "mass " << ~(piminus + piplus) << " " << ~isobar2 << " " << isobar2_mass << " " << piplus_p << " " << piplus.r() << std::endl;
	  std::cout << std::endl << std::endl;
	}



    
    // pi0 decay to 2 gamma, in pi0 rest frame
    double gam1_p = CMmomentum (PI0_MASS, 0.0, 0.0);
    gamma1.polar(gam1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    gamma1.t(gam1_p);
    gamma2.set(gam1_p,threeVec (0.0, 0.0, 0.0) - gamma1.V());
    gamma2.t(gam1_p);
    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar2_p * piplus_p * phi_p * Kp_p * gam1_p );
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;
        
         // Get K's back
    
    tmp.set(phi.t(),zeroVec - phi.V());
        Boost.set (tmp);
        Kp = Boost * Kp;
	Km = Boost * Km; 

	// boost gammas to omega  rest frame
	tmp.set(pizero.t(),zeroVec - pizero.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;

	// boost pi+ pi- to omega rest frame
	tmp.set(isobar2.t(),zeroVec - isobar2.V());
	Boost.set(tmp);



	piplus = Boost * piplus;
	piminus = Boost * piminus;

	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi+ (isobar2) : "; pParticleGamp(PiPlus,piplus);
	  std::cout << "pi- (isobar2) : "; pParticleGamp(PiMinus,piminus);
	  std::cout << "mass " << ~(piminus + piplus) << " " << ~isobar2 << " " << isobar2_mass << std::endl;
	  std::cout << std::endl << std::endl;
	}



	// boost omega progeny
	tmp.set(omeg.t(),zeroVec - omeg.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	piplus = Boost * piplus;
	piminus = Boost * piminus;
	pizero = Boost * pizero;
	isobar2 = Boost * isobar2;


	

	if (debug) {
	  std::cout << "\n";
	  std::cout << "Kp : "; pParticleGamp(KPlus,Kp);
	  std::cout << "Km : "; pParticleGamp(KMinus,Km);
	  std::cout << std::endl << std::endl;
	}

        // boost  everyone to resonance   frame

        tmp.set(resonance.t(),zeroVec - resonance.V());
        Boost.set (tmp);
        piplus = Boost * piplus;
	piminus = Boost * piminus;
        pizero = Boost * pizero;
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	Kp = Boost * Kp;
	Km = Boost * Km;
        phi = Boost * phi;
	isobar2 = Boost * isobar2;
	omeg = Boost * omeg;


        // boost from CM to target rest frame (lab)
        tmp.set (CMtarget.t(),zeroVec - CMtarget.V());
	Boost.set(CMtarget);
	CMtarget = Boost * CMtarget;
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        piplus = Boost * piplus;
	piminus = Boost * piminus;
        pizero = Boost * pizero;
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	Kp = Boost * Kp;
	Km = Boost * Km;
        phi = Boost * phi;
	omeg = Boost * omeg;
	isobar2 = Boost * isobar2;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "pi+ :\n  ";
	  piplus.print();
	  std::cerr << "pi-2:\n  ";
	  pizero.print();
	  std::cerr << "pi-1:\n  ";
	  piminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << ~(pizero + piminus + piplus) << " ";
	  std::cout << pow (~(pizero + piplus), 2.0) << " "
	       << pow (~(pizero + piminus), 2.0) << " "
	       << pow (~(piplus + piminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << (pizero + piminus + piplus + recoil).lenSq() << " ";
	  std::cout << std::endl;
	} 
	else if (txt2part_style){
	  if (printGamma) {
	    std::cout << "7" << std::endl;
	  } 
	  else {
	    std::cout << "6" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Baryon, production, recoil);
	  pParticle_txt2part(KPlus, production, Kp);
	  pParticle_txt2part(KMinus, production, Km);
	  pParticle_txt2part(PiMinus,production,piminus);
	  pParticle_txt2part(PiPlus,production,piplus);
	  if (printGamma){
	    pParticle_txt2part(Gamma, production, gamma1);
	    pParticle_txt2part(Gamma, production, gamma2);
	  }
	  else{
	    pParticle_txt2part(Pi0,production,pizero);
	  }
	}
	else if (gamp){
	  if (printGamma) {
	    std::cout << "8" << std::endl;
	  } 
	  else {
	    std::cout << "7" << std::endl;
	  }
	  pParticleGamp(Beam,beam);
	  pParticleGamp(Baryon,recoil);
	  pParticleGamp(KPlus, Kp);
	  pParticleGamp(KMinus,Km);
	  pParticleGamp(PiMinus,piminus);
	  pParticleGamp(PiPlus,piplus);
	  if (printGamma){
	    pParticleGamp(Gamma, gamma1);
	    pParticleGamp(Gamma, gamma2);
	  }
	  else{
	    pParticleGamp(Pi0,pizero);
	  }
	}
	else if (debug) {
	  std::cout << "5" << std::endl;
	  std::cout << "beam "; pParticleGamp(Beam,beam);
	  std::cout << "recoil "; pParticleGamp(Baryon,recoil);
	  std::cout << "Resonance "; pParticleGamp(Rho0,resonance);
	  std::cout << "omega "; pParticleGamp(omega,omeg);
	  std::cout << "isobar2 "; pParticleGamp(Rho0,isobar2);
	  std::cout << "pizero (omega)"; pParticleGamp(Pi0,pizero);
	  std::cout << "piplus (omega) "; pParticleGamp(PiPlus,piplus);
	  std::cout << "piminus (omega) "; pParticleGamp(PiMinus,piminus);
	}
	  
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Beam,production,beam);


	  if (printBaryon)
	    pParticle(Baryon,production,recoil);
	  pParticle(KPlus,production,Kp);
	  pParticle(KMinus,production,Km);
	  pParticle(PiPlus,production,piminus);
	  pParticle(PiMinus,production,piplus);
	  if (printGamma) {
	    pParticle(Gamma,production,gamma1);
	    pParticle(Gamma,production,gamma2);
	  }
	  else {
	    pParticle(Pi0,production,pizero);
	  }
	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}







/* double Dalitz */
void doubleDalitz (int argc, char *argv[],Particle_t Beam,Particle_t Pi1,Particle_t Pi2,Particle_t Recoil,int decay)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  Particle_t Baryon;
  double 
    masslow = Mass(Pi1) + Mass(Pi2),
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    pi1,
    pi2,
    gamma1,
    gamma2,
    em1,em2,ep1,ep2,
    isobar1;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass;
  float baryonmass;
  int printBaryon = 0;
  int printGamma = 0;
  int debug = 0;
  Particle_t Target = Proton;
  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break; 
      case 'D':
	debug = 1;
	break;
      case 'h':
	UsagedoubleDalitz(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }
 
  double slope = (Slope > 0.0) ? Slope : 10.0;
  double pipiThreshold = Mass(Pi1) + Mass(Pi2);
  int qtot = Q(Beam) + Q(Proton) - Q(Pi1) - Q(Pi2);
  switch (qtot) {
  case 0:
    Baryon = Neutron;
    break;
  case 1:
    Baryon = Proton;
    break;
  default:
    std::cerr << "illegal charge combination:\ttotal charge  = " << qtot << std::endl;
    exit(0);
  }


  beamMass = Mass(Beam);
  baryonmass = Mass(Baryon);

  if(masslow < pipiThreshold)
    masslow = pipiThreshold;

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - Mass(Baryon);
    if (masshigh < pipiThreshold) {
      std::cerr << "Meson high mass below 2 particle threshold" << std::endl;
      exit(1);
    }
    else if (masshigh > ( CMenergy - Mass(Baryon))) {
      masshigh =  CMenergy - Mass(Baryon);
    }

    do {
      resonance_mass = randm(masslow*masslow, masshigh*masshigh);
      resonance_mass = sqrt(resonance_mass);
    } while (resonance_mass < pipiThreshold);
    
    double e1_p,e2_p,gamma1_p,gamma2_p,pi_p;
    double beam_p      = CMmomentum (CMenergy, Mass(Beam), PROTON_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, Mass(Baryon));
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,Mass(Beam),PROTON_MASS,resonance_mass,Mass(Baryon));
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      twoBodyDecay(CMenergy,resonance_mass,Mass(Baryon),resonance,recoil);

    else {
      do {
	t = getT(tMin,slope);
	costheta = 1 + (t - tMin)/(2 * resonance_p * beam_p);
	if (debug) {
	  cout << "T " << t << " " << tMin << " " << slope << " " << costheta << endl;
	}
      } while (fabs(costheta) >=  1.0);

    
      resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
      resonance.t(resonance_E);
      recoil.set(sqrt(resonance.V().lenSq() +  pow (Mass(Baryon), 2.0)),zeroVec - resonance.V());
    }


    /*
     *  now do pi0 pi0 decay in resonance rest frame
     */

    twoBodyDecay(resonance_mass,Mass(Pi0),Mass(Pi0),pi1,pi2);
    pi_p = CMmomentum(resonance_mass,Mass(Pi1),Mass(Pi2));

    // neutral decay to gamma e+ e-, in pi0 rest frame
#define EPEMSLOPE 70.0
    if (decay) {
      fourVec tmp;

      fourVec x;
      double epemMass;
      do {
	-getT(-2 * Mass(Electron),EPEMSLOPE);
      } while (epemMass > Mass(Pi1));
      twoBodyDecay(Mass(Pi1),0.0,epemMass,gamma1,x);
      gamma1_p = CMmomentum(Mass(Pi1),epemMass,Mass(Gamma));
     
      if (debug) {
	cout << "XXX " << ~(gamma1 + x) << " " << ~gamma1 << endl;
      }
      // x and gamma1 are in pi01 rest frame
      twoBodyDecay(epemMass,Mass(Electron),Mass(Positron),em1,ep1); 
      e1_p = em2.V().r();
      if (debug) {
	cout << "DBG " << epemMass << " " << ~(ep1 + em1) << " ";
      }
      // ep1 and em1 are in x rest frame, x is in rest frame of pi01
      tmp.set(x.t(),zeroVec - x.V());
      Boost.set (tmp);
      ep1 *= Boost;
      em1 *= Boost;
      // ep1 and em1 are now in frame of the parent pi01, as is x and gamma1
      if (debug) {
	fourVec X = x - ep1 - em1;
	cout << ~x  << " " << ~(ep1 + em1) << " ";
	cout << X.V().r() << " " << Mass(Pi1) << " ";
	cout << ~(gamma1 + x) << " ";
      }
      tmp.set(pi1.t(),zeroVec - pi1.V());
      Boost.set(tmp);
      ep1 *= Boost;
      em1 *= Boost;
      x *= Boost;
      gamma1 *= Boost;

      // now everything (x,gamma1,ep1,em1)  is in pi0 pi0  resonance rest frame

      if (debug) {
	fourVec Xee = ep1 + em1;
	fourVec X = x - Xee;
	fourVec XXX = gamma1 + Xee;
	fourVec X2 = pi1 - XXX;
	cout << X.V().r() << " " << X2.V().r() << " " << " " << ~(gamma1 + x) << " " << ~(em1 + ep1) << " ";
      }
      
      do {
	epemMass = -getT(-2 * Mass(Electron),EPEMSLOPE);
      } while(epemMass > Mass(Pi2));
      twoBodyDecay(Mass(Pi2),0.0,epemMass,gamma2,x);
      twoBodyDecay(epemMass,Mass(Electron),Mass(Positron),em2,ep2); 
      gamma2_p = CMmomentum(Mass(Pi2),epemMass,Mass(Gamma));
      e2_p = em2.V().r();
      tmp.set(x.t(),zeroVec - x.V());
      Boost.set (tmp);
      ep2 *= Boost;
      em2 *= Boost;  
      tmp.set(pi2.t(),zeroVec - pi2.V());
      Boost.set(tmp);
      ep2 *= Boost;
      em2 *= Boost;
      x *= Boost;
      gamma2 *= Boost;


      if (debug) {
	fourVec Xee = ep2 + em2;
	fourVec X = x - Xee;
	fourVec XXX = gamma2 + Xee;
	fourVec X2 = pi2 - XXX;
	cout << epemMass << " " << ~(gamma2 + x) << " " << ~(em2 + ep2) << " ";
	cout << X.V().r() << " " << X2.V().r() << " " << " " << ~(gamma2 + x) << " " << ~(em2 + ep2) << " ";
      }
      // first transform the pi0s
      tmp.set(resonance.t(),zeroVec - resonance.V());
      Boost.set (tmp);
      pi1 *= Boost;
      pi2 *= Boost; 
      ep1 *= Boost;
      em1 *= Boost;
      ep2 *= Boost;
      em2 *= Boost;
      gamma1 *= Boost;
      gamma2 *= Boost;
	
      // boost from CM to target rest frame (lab)
      Boost.set (CMtarget);
      resonance *= Boost;
      recoil *= Boost;
      pi1 *= Boost;
      pi2 *= Boost; 
      ep1 *= Boost;
      em1 *= Boost;
      ep2 *= Boost;
      em2 *= Boost;
      gamma1 *= Boost;
      gamma2 *= Boost;
    }

    if (debug) {
      cout << ~(em1 + ep1) << " " <<  ~(em2 + ep2) << " " << ~(em1 + ep1)/PI0_MASS << " " << ~(ep2 + em2)/PI0_MASS << endl;
    }


   /*
     *  compute lorentz factor
     */
    LorentzFactor = resonance_p * pi_p * gamma1_p * gamma2_p * e1_p * e2_p;
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {



    // write out


    if (gamp) {

      event evt;
	  
      particle gBeam(PDGtable.get(pid2name(Beam)),Q(Beam));
      particle gTarget(PDGtable.get(pid2name(Target)),Q(Target));
      particle gRecoil(PDGtable.get(pid2name(Baryon)),Q(Baryon));
      particle gPi1(PDGtable.get(pid2name(Pi1)),Q(Pi1));
      particle gPi2(PDGtable.get(pid2name(Pi2)),Q(Pi2));
      particle gGamma1(PDGtable.get(pid2name(Gamma)),Q(Gamma));
      particle gGamma2(PDGtable.get(pid2name(Gamma)),Q(Gamma));

      particle gEp1(PDGtable.get(pid2name(Positron)),Q(Positron));
      particle gEm1(PDGtable.get(pid2name(Electron)),Q(Electron));
      particle gEp2(PDGtable.get(pid2name(Positron)),Q(Positron));
      particle gEm2(PDGtable.get(pid2name(Electron)),Q(Electron));


      gBeam.set4P(beam);
      gTarget.set4P(target);
      gRecoil.set4P(recoil);
      gPi1.set4P(pi1);
      gPi2.set4P(pi2);
      if (decay) {
	gGamma1.set4P(gamma1);
	gGamma2.set4P(gamma2);
	gEp1.set4P(ep1);
	gEm1.set4P(em1);
	gEp2.set4P(ep2);
	gEm2.set4P(em2);
      }

      evt.beam(gBeam);
      evt.target(gTarget);
      evt.addfinal(gRecoil);
      evt.addfinal(gPi1);
      evt.addfinal(gPi2);
      if (decay) {
	evt.addfinal(gGamma1);
	evt.addfinal(gGamma2);
	evt.addfinal(gEp1);
	evt.addfinal(gEm1);
	evt.addfinal(gEp2);
	evt.addfinal(gEm2);
      }


      std::cout << evt;
    }
    else if (txt2part_style){
#define PIPIDALITZ 22
      txtEvent evt; 
      particle gBeam(PDGtable.get(pid2name(Beam)),Q(Beam));
      particle gTarget(PDGtable.get(pid2name(Target)),Q(Target));
      particle gRecoil(PDGtable.get(pid2name(Baryon)),Q(Baryon));
      particle gPi1(PDGtable.get(pid2name(Pi1)),Q(Pi1));
      particle gPi2(PDGtable.get(pid2name(Pi2)),Q(Pi2));
      particle gGamma1(PDGtable.get(pid2name(Gamma)),Q(Gamma));
      particle gGamma2(PDGtable.get(pid2name(Gamma)),Q(Gamma));

      particle gEp1(PDGtable.get(pid2name(Positron)),Q(Positron));
      particle gEm1(PDGtable.get(pid2name(Electron)),Q(Electron));
      particle gEp2(PDGtable.get(pid2name(Positron)),Q(Positron));
      particle gEm2(PDGtable.get(pid2name(Electron)),Q(Electron));


      gBeam.set4P(beam);
      gTarget.set4P(target);
      gRecoil.set4P(recoil);
      gPi1.set4P(pi1);
      gPi2.set4P(pi2);
      if (decay) {
	gGamma1.set4P(gamma1);
	gGamma2.set4P(gamma2);
	gEp1.set4P(ep1);
	gEm1.set4P(em1);
	gEp2.set4P(ep2);
	gEm2.set4P(em2);
      }

      evt.beam(gBeam);
      evt.target(gTarget);
      evt.addfinal(gRecoil);
      evt.addfinal(gPi1);
      evt.addfinal(gPi2);
      if (decay) {
	evt.addfinal(gGamma1);
	evt.addfinal(gGamma2);
	evt.addfinal(gEp1);
	evt.addfinal(gEm1);
	evt.addfinal(gEp2);
	evt.addfinal(gEm2);
      }
      evt.code(PIPIDALITZ);
      evt.weight(1.0);

      std::cout << evt;
    }
    maxevents--;

      }
    }
  }
}

void twoBodyDecay(double mass,double mass1,double mass2,fourVec &p1,fourVec &p2)
{
  double cmp =  CMmomentum (mass,mass1, mass2);
  p1.polar(cmp,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
  p1.t(sqrt(mass1*mass1 + cmp*cmp));
  p2.set(cmp,threeVec (0.0, 0.0, 0.0) - p1.V());
  p2.t(sqrt(cmp*cmp + mass2*mass2));
 
}
	   	   



/* ---------------- beam p -> c+ c- n0 N -------------------------------------*/

/*  example pi- p -> K+ K- eta n */

void cpcmn0N (int argc, char *argv[],Particle_t Beam,Particle_t Cplus,Particle_t Cminus,Particle_t N0,int decay)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  int PrintGammas = 0;
  Particle_t Baryon;
  double 
    masslow = Mass(Cplus) + Mass(Cminus) + Mass(N0),
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    cminus,
    cplus,
    n0,
    gamma1,
    gamma2,
    isobar1;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass;
  float baryonmass;
  int printBaryon = 0;
  int printGamma = 0;
  int debug = 0;
  Particle_t Target = Proton;
  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break; 
      case 'D':
	debug = 1;
	break;
      case 'h':
	Usagecpcmn0N(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 10.0;
  double cpcmn0Threshold = Mass(Cplus) + Mass(Cminus) + Mass(N0);
  double cpcmThreshold = Mass(Cplus) + Mass(Cminus); 
  int qtot = Q(Beam) + Q(Proton) - Q(Cplus) - Q(Cminus) -Q(N0);
  switch (qtot) {
  case 0:
    Baryon = Neutron;
    break;
  case 1:
    Baryon = Proton;
    break;
  default:
    std::cerr << "illegal charge combination:\ttotal charge  = " << qtot << std::endl;
    exit(0);
  }


  beamMass = Mass(Beam);
  baryonmass = Mass(Baryon);

  if(masslow < cpcmn0Threshold)
    masslow = cpcmn0Threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - Mass(Baryon);
    if (masshigh < cpcmn0Threshold) {
      std::cerr << "Meson high mass below 3 particle threshold" << std::endl;
      exit(1);
    }
    else if (masshigh > ( CMenergy - Mass(Baryon))) {
      masshigh =  CMenergy - Mass(Baryon);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < cpcmn0Threshold);
    isobar1_mass = randm(cpcmThreshold, (resonance_mass - Mass(N0)));

    double beam_p      = CMmomentum (CMenergy, Mass(Beam), PROTON_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, Mass(Baryon));
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,Mass(Beam),PROTON_MASS,resonance_mass,Mass(Baryon));
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);

    else {
        costheta = cosThetaExp(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass,slope);
    //     costheta = 1.0;
    }
    

    if ((resonance_p < 0.0) || (costheta > 1) || (costheta < -1))  {
      std::cerr << "Error" << std::endl;
    }

    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (Mass(Baryon), 2.0)),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, Mass(N0));
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

    isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    n0.set(sqrt(isobar1.V().lenSq() + pow (Mass(N0), 2.0)),zeroVec -isobar1.V());




    /*
     *  now do decay in isobar1 rest frame
     */
    // c+ c-
    double cplus_p = CMmomentum (isobar1_mass, Mass(Cplus),Mass(Cminus));

    cplus.polar( cplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   cplus.polar( cplus_p,acos (1.0),0.0);
    cplus.t(sqrt (cplus.V().lenSq () + pow (Mass(Cplus), 2.0)));
    cminus.set(sqrt (cplus.V().lenSq () + pow (Mass(Cminus), 2.0)),zeroVec - cplus.V());
    
    // neutral decay to 2gamma, in pi0 rest frame
    if (decay) {
    double gam1_p = CMmomentum (Mass(N0), 0.0, 0.0);
    gamma1.polar(gam1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    gamma1.t(gam1_p);
    gamma2.set(gam1_p,threeVec (0.0, 0.0, 0.0) - gamma1.V());
    gamma2.t(gam1_p);
    }
    /*
     *  compute lorentz factor
     */
    LorentzFactor = FlatMass ? 1 : (isobar1_p * cplus_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
	  /* generate vertices */
	  threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

	// boost gammas to neutral rest frame
	if (decay) {
	tmp.set(n0.t(),zeroVec - n0.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	}
        // boost from c+ c- rest frame to c+ c- n0(X) rest frame

        tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        cplus = Boost * cplus;
        cminus = Boost * cminus;

        // boost from 3 particle meson  rest frame to CM
        tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);

        isobar1 = Boost * isobar1;
        cplus = Boost * cplus;
	cminus = Boost * cminus;
        n0 = Boost * n0;
	if (decay) {
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	}
 

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        isobar1 = Boost * isobar1;
        cplus = Boost * cplus;
        cminus = Boost * cminus;
        n0 = Boost * n0;
	if (decay) {
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	}



	
	

	if (Print) {

	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "isobar1\n";
	  std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	  std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	  std::cerr << "  isobar1:\n    ";
	  isobar1.print();
	  std::cerr << "c+ :\n  ";
	  cplus.print();
	  std::cerr << "n0:\n  ";
	  n0.print();
	  std::cerr << "c-:\n  ";
	  cminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << resonance_mass << " " << isobar1_mass*isobar1_mass  << " ";
	  std::cout << ~(n0 + cminus + cplus) << " ";
	  std::cout << pow (~(n0 + cplus), 2.0) << " "
	       << pow (~(n0 + cminus), 2.0) << " "
	       << pow (~(cplus + cminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << (n0 + cminus + cplus + recoil).lenSq() << " ";
	  std::cout << (beam - cminus - cplus - gamma1 - gamma2).lenSq() << " ";
	  std::cout << std::endl;
	} 
	else if (debug) {
	  std::cout << "Z " << production.x() << " " << production.y() << " " << production.z() << std::endl;
	
	}
	else if (txt2part_style){
	  if (decay && printGamma) {
	    std::cout << "5" << std::endl;
	  } 
	  else {
	    std::cout << "4" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Baryon, production, recoil);
	  pParticle_txt2part(Cminus, production, cminus);
	  pParticle_txt2part(Cplus, production, cplus);
	  if (decay && printGamma){
	    pParticle_txt2part(Gamma, production, gamma1);
	    pParticle_txt2part(Gamma, production, gamma2);
	  }
	  else{
	    pParticle_txt2part(N0,production,n0);
	  }
	}
	else if (gamp) {

	event evt;
	  
	particle gBeam(PDGtable.get(pid2name(Beam)),Q(Beam));
	particle gTarget(PDGtable.get(pid2name(Target)),Q(Target));
	particle gRecoil(PDGtable.get(pid2name(Baryon)),Q(Baryon));
	particle gCminus(PDGtable.get(pid2name(Cminus)),Q(Cminus));
	particle gCplus(PDGtable.get(pid2name(Cplus)),Q(Cplus));
	particle gGamma1(PDGtable.get(pid2name(Gamma)),Q(Gamma));
	particle gGamma2(PDGtable.get(pid2name(Gamma)),Q(Gamma));
	particle gN0(PDGtable.get(pid2name(N0)),Q(N0));

	gBeam.set4P(beam);
	gTarget.set4P(target);
	gRecoil.set4P(recoil);
	gCminus.set4P(cminus);
	gCplus.set4P(cplus);
	if (decay) {
	gGamma1.set4P(gamma1);
	gGamma2.set4P(gamma2);
	}
	gN0.set4P(n0);

	evt.beam(gBeam);
	evt.target(gTarget);
	evt.addfinal(gRecoil);
	evt.addfinal(gCminus);
	evt.addfinal(gCplus);
	if (decay && printGamma) {
	  evt.addfinal(gGamma1);
	  evt.addfinal(gGamma2);
	}
	else if (printAll) {
	  if (decay) {
	  evt.addfinal(gGamma1);
	  evt.addfinal(gGamma2);
	  }
	  evt.addfinal(gN0);
	}
	else {
	  evt.addfinal(gN0);
	}
	  
	std::cout << evt;

	}
      else {	  	  
	/*
	 *  write event
	 */
	

	if (printBeam)
	  pParticle(Beam,production,beam);


	if (printBaryon)
	  pParticle(Proton,production,recoil);
	pParticle(Cplus,production,cminus);
	pParticle(Cminus,production,cplus);
	if (decay && printGamma) {
	  pParticle(Gamma,production,gamma1);
	  pParticle(Gamma,production,gamma2);
	}
	else {
	  pParticle(N0,production,n0);
	}
      }
      nevents++;
      if (verbose) {
	if (!(nevents % 100)) 
	  std::cerr << nevents << "\r" << flush;
      }
      maxevents--;
      }

    
    }
  }

  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

/*----------------End of charged1 charged2 neutral-----------------------------------*/


/* ---------------- beam p -> c+ c- n0 N -------------------------------------*/

/*  example pi- p -> K+ K- n */

void cpcm(int argc, char *argv[],Particle_t Beam,Particle_t Cplus,Particle_t Cminus,Particle_t Recoil)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  Particle_t Baryon;
  double 
    masslow = Mass(Cplus) + Mass(Cminus),
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    cminus,
    cplus;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass;
  float baryonmass;
  int printBaryon = 0;
  int debug = 0;  

  Particle_t Target = Proton;

  /* generate vertices */
  threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'D':
	debug = 1;
	break;
      case 'h':
	Usagecpcm(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 10.0;
  int qtot = Q(Beam) + Q(Proton) - Q(Cplus) - Q(Cminus);
  double cpcmThreshold =  Mass(Cplus) + Mass(Cminus);

  Baryon = Recoil;

  beamMass = Mass(Beam);
  baryonmass = Mass(Baryon);

  if(masslow < cpcmThreshold)
    masslow = cpcmThreshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - Mass(Baryon);
    if (masshigh < cpcmThreshold) {
      std::cerr << "Meson high mass below 2 particle threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < cpcmThreshold);

    double beam_p      = CMmomentum (CMenergy, Mass(Beam), PROTON_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, Mass(Baryon));
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double costhetax,tx;
    double t;
    double mTot = beamMass + TARGET_MASS + Mass(Baryon) + Mass(Cplus) - Mass(Cminus);
    double s =  pow(sqrt(beam_p * beam_p + beamMass * beamMass)  + sqrt(beam_p * beam_p + TARGET_MASS *TARGET_MASS),2);

    tMin = tmin(pbeam.z,Mass(Beam),PROTON_MASS,resonance_mass,Mass(Baryon));

    
    
    /* use distribution t' * exp(-slope*t'), where t' = abs(t - tmin) and
       0 <= t' <= 4*pbeam*presonance in CM */
    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass);
    //    std::cerr << tMin << std::endl;
    t_max = 4. * beam_p * resonance_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {

        double tp = tprimeExp(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass,slope);
        costheta = 1. - 2.*t/t_max;
        costheta = cosThetaExp(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass,slope);
    //     costheta = 1.0;
    }




    if (uChannel)
      resonance.polar(resonance_p,M_PI + acos (costheta), randm (-M_PI,M_PI));
    else
      resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (Mass(Baryon), 2.0)),zeroVec - resonance.V());


    /*
     *  now do decay in isobar1 rest frame
     */
    // c+ c-
    double cplus_p = CMmomentum (resonance_mass, Mass(Cplus),Mass(Cminus));

    cplus.polar( cplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //     cplus.polar( cplus_p,acos (1.0),0.0);
    cplus.t(sqrt (cplus.V().lenSq () + pow (Mass(Cplus), 2.0)));
    cminus.set(sqrt (cplus.V().lenSq () + pow (Mass(Cminus), 2.0)),zeroVec - cplus.V());
    

    /*
     *  compute lorentz factor
     */
    LorentzFactor = FlatMass ? 1 : resonance_p * cplus_p;
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;


        // boost from c+ c- rest frame to CM rest frame

        tmp.set(resonance.t(),zeroVec - resonance.V());
        Boost.set (tmp);
        cplus = Boost * cplus;
        cminus = Boost * cminus;


        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        cplus = Boost * cplus;
        cminus = Boost * cminus;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "c+ :\n  ";
	  cplus.print();
	  std::cerr << "c-:\n  ";
	  cminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << "DALITZ ";
	  std::cout << resonance_mass << " ";
	  std::cout << ~(recoil + cminus + cplus) << " ";
	  std::cout << pow (~(recoil + cplus), 2.0) << " "
	       << pow (~(recoil + cminus), 2.0) << " "
	       << pow (~(cplus + cminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << std::endl;
	} 
	else if (debug) {
	  std::cout << "T " << (beam - cplus - cminus).lenSq() << " " << (target - recoil).lenSq() << " " << t << " " << tMin << " " << costheta << " " <<  costhetax << " " << (t - tMin) << " " << tx << " " << std::endl;
	}
	if (txt2part_style){
	  std::cout << "3" << std::endl;
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Recoil, production, recoil);
	  pParticle_txt2part(Cminus, production, cminus);
	  pParticle_txt2part(Cplus, production, cplus);
	}
	else if (gamp) {
	  event evt;
	  
	  particle gBeam(PDGtable.get(pid2name(Beam)),Q(Beam));
	  particle gTarget(PDGtable.get(pid2name(Target)),Q(Target));
	  particle gRecoil(PDGtable.get(pid2name(Recoil)),Q(Recoil));
	  particle gCminus(PDGtable.get(pid2name(Cminus)),Q(Cminus));
	  particle gCplus(PDGtable.get(pid2name(Cplus)),Q(Cplus));

	  gBeam.set4P(beam);
	  gTarget.set4P(target);
	  gRecoil.set4P(recoil);
	  gCminus.set4P(cminus);
	  gCplus.set4P(cplus);

	  evt.beam(gBeam);
	  evt.target(gTarget);
	  evt.addfinal(gRecoil);
	  evt.addfinal(gCminus);
	  evt.addfinal(gCplus);
	  
	  std::cout << evt;
	  
	  
	}
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Beam,production,beam);


	  pParticle(Recoil,production,recoil);
	  pParticle(Cplus,production,cminus);
	  pParticle(Cminus,production,cplus);

	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

/*----------------End of charged1 charged2 -----------------------------------*/







/*  etaprime pi */

void etaprimepi(int argc, char *argv[],Particle_t Beam,Particle_t Pi)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  int PrintGammas = 0;
  Particle_t Baryon;
  double 
    masslow = ETAPRIME_MASS + PI_MASS,
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    etaprime_p;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    etaprime,
    pi,
    piplus,
    piminus,
    eta,
    gamma1,
    gamma2;
  fourVec isobar1,isobar2;
  // isobar1 is pi pi isobar off of omega pi pi system
  // isobar2 is pi+ pi- isobar from omega decay
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = GAMMA_MASS;
  int printBaryon = 0;
  int printGamma = 0;

  double pimass  = Mass(Pi),baryonmass = Mass(Baryon);

  float tMin;

  int debug = 0;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'D':
	debug = 1;
	break;
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break;
      case 'h':
	UsageM26(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 2.9;
  double threshold = ETAPRIME_MASS + pimass;

  int qtot = Q(Beam) + Q(Proton) - Q(Pi);

  switch (qtot) {
  case 0:
    Baryon = Neutron;
    break;
  case 1:
    Baryon = Proton;
    break;
  default:
    std::cerr << "illegal charge combination:\ttotal charge  = " << qtot << std::endl;
    exit(0);
}


  beamMass = Mass(Beam);
  baryonmass = Mass(Baryon);

  // set up the recoil

  
  if(masslow < threshold)
    masslow = threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - baryonmass;
    if (masshigh < threshold) {
      std::cerr << "resonance  high mass below threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < threshold);

    double beam_p      = CMmomentum (CMenergy, beamMass, TARGET_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, baryonmass);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    tMin = tmin(pbeam.z,beamMass,PROTON_MASS,resonance_mass,baryonmass);


    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {
      t = getT(tMin,slope);
      costheta = 1 + (t - tMin)/(2 * resonance_p * beam_p);
    }

  
    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (baryonmass, 2.0)),zeroVec - resonance.V());



    /*
     *  now do decay in resonance rest frame
     */
     
    etaprime_p = CMmomentum(resonance_mass,ETAPRIME_MASS,pimass);
    etaprime.polar(etaprime_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI));
    etaprime.t(sqrt (etaprime.V().lenSq () + pow (ETAPRIME_MASS, 2.0)));
    pi.set(sqrt (etaprime.V().lenSq () + pow (PI_MASS, 2.0)),zeroVec - etaprime.V());

     // code to do etaprime decay
    double himass = ETAPRIME_MASS - ETA_MASS;
    double lomass = 2 * PI_MASS;
    double isobar2_mass = randm(lomass,himass);
    double isobar2_p = CMmomentum(ETAPRIME_MASS,isobar2_mass,ETA_MASS);
    isobar2.polar( isobar2_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI));
    isobar2.t(sqrt (isobar2.V().lenSq () + pow (isobar2_mass, 2.0))); 
    eta.set(sqrt(isobar2.V().lenSq() + pow(ETA_MASS,2.0)),zeroVec - isobar2.V());

 
    // isobar2 decay
    double piplus_p = CMmomentum(isobar2_mass,PI_MASS,PI_MASS);

    piplus.polar( piplus_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    piplus.t(sqrt(piplus.V().lenSq() + pow(PI_MASS,2.0)));
    piminus.set(sqrt(piplus.V().lenSq() + pow(PI_MASS,2.0)),zeroVec - piplus.V());
	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi+ (isobar2) : "; pParticleGamp(PiPlus,piplus);
	  std::cout << "pi- (isobar2) : "; pParticleGamp(PiMinus,piminus);
	  std::cout << "mass " << ~(piminus + piplus) << " " << ~isobar2 << " " << isobar2_mass << " " << piplus_p << " " << piplus.r() << std::endl;
	  std::cout << std::endl << std::endl;
	}



    
    // eta decay to 2 gamma, in eta rest frame
    double gam1_p = CMmomentum (ETA_MASS, 0.0, 0.0);
    gamma1.polar(gam1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    gamma1.t(gam1_p);
    gamma2.set(gam1_p,threeVec (0.0, 0.0, 0.0) - gamma1.V());
    gamma2.t(gam1_p);
    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar2_p * etaprime_p * piplus_p * gam1_p );
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

	// boost gammas to eta rest frame
	tmp.set(eta.t(),zeroVec - eta.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;

	// boost pi+ pi-  to isobar2  rest frame
	tmp.set(isobar2.t(),zeroVec - isobar2.V());
	Boost.set(tmp);

	piplus = Boost * piplus;
	piminus = Boost * piminus;

	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi+ (isobar2) : "; pParticleGamp(PiPlus,piplus);
	  std::cout << "pi- (isobar2) : "; pParticleGamp(PiMinus,piminus);
	  std::cout << "mass " << ~(piminus + piplus) << " " << ~isobar2 << " " << isobar2_mass << std::endl;
	  std::cout << std::endl << std::endl;
	}



	// boost etaprime  progeny
	tmp.set(etaprime.t(),zeroVec - etaprime.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	piplus = Boost * piplus;
	piminus = Boost * piminus;
	eta = Boost * eta;
	isobar2 = Boost * isobar2;
	if (debug) {
	  std::cout << "\n";
	  std::cout << "pi : "; pParticleGamp(Pi,pi);
	  std::cout << std::endl << std::endl;
	}

        // boost  everyone to resonance   frame

        tmp.set(resonance.t(),zeroVec - resonance.V());
        Boost.set (tmp);
        piplus = Boost * piplus;
	piminus = Boost * piminus;
        eta = Boost * eta;
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	pi = Boost * pi;
	isobar2 = Boost * isobar2;
	etaprime = Boost * etaprime;

        // boost from CM to target rest frame (lab)
        tmp.set (CMtarget.t(),zeroVec - CMtarget.V());
	Boost.set(CMtarget);
	CMtarget = Boost * CMtarget;
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        piplus = Boost * piplus;
	piminus = Boost * piminus;
        eta = Boost * eta;
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	pi = Boost * pi;
	etaprime = Boost * etaprime;
	isobar2 = Boost * isobar2;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "pi+ :\n  ";
	  piplus.print();
	  std::cerr << "eta:\n  ";
	  eta.print();
	  std::cerr << "pi- ((from etprime):\n  ";
	  piminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << resonance_mass << " ";
	  std::cout << ~(eta + piminus + piplus) << " ";
	  std::cout << pow (~(eta + piplus), 2.0) << " "
	       << pow (~(eta + piminus), 2.0) << " "
	       << pow (~(piplus + piminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << (eta + piminus + piplus + recoil).lenSq() << " ";
	  std::cout << std::endl;
	} 
	else if (txt2part_style){
	  if (printGamma) {
	    std::cout << "7" << std::endl;
	  } 
	  else {
	    std::cout << "6" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Baryon, production, recoil);
	  pParticle_txt2part(Pi, production, pi);
	  pParticle_txt2part(PiMinus,production,piminus);
	  pParticle_txt2part(PiPlus,production,piplus);
	  if (printGamma){
	    pParticle_txt2part(Gamma, production, gamma1);
	    pParticle_txt2part(Gamma, production, gamma2);
	  }
	  else{
	    pParticle_txt2part(Eta,production,eta);
	  }
	}
	else if (gamp){
	  if (printGamma || printAll) {
	    std::cout << "9" << std::endl;
	  } 
	  else {
	    std::cout << "7" << std::endl;
	  }
	  pParticleGamp(Beam,beam);
	  pParticleGamp(Baryon,recoil);
	  pParticleGamp(Pi, pi);
	  pParticleGamp(PiMinus,piminus);
	  pParticleGamp(PiPlus,piplus);
	  pParticleGamp(EtaPrime,etaprime);
	  pParticleGamp(Eta,eta);
	  if (printGamma || printAll){
	    pParticleGamp(Gamma, gamma1);
	    pParticleGamp(Gamma, gamma2);
	  }
	}
	else if (debug) {
	  std::cout << "5" << std::endl;
	  std::cout << "beam "; pParticleGamp(Beam,beam);
	  std::cout << "recoil "; pParticleGamp(Baryon,recoil);
	  std::cout << "Resonance "; pParticleGamp(Rho0,resonance);
	  std::cout << "etaprime "; pParticleGamp(EtaPrime,etaprime);
	  std::cout << "isobar2 "; pParticleGamp(Rho0,isobar2);
	  std::cout << "pi (isobar1) "; pParticleGamp(Pi,pi);
	  std::cout << "eta (etaprime)"; pParticleGamp(Eta,eta);
	  std::cout << "piplus (etaprime) "; pParticleGamp(PiPlus,piplus);
	  std::cout << "piminus (etaprime) "; pParticleGamp(PiMinus,piminus);
	}
	  
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Beam,production,beam);


	  if (printBaryon)
	    pParticle(Baryon,production,recoil);
	  pParticle(Pi,production,pi);
	  pParticle(PiPlus,production,piminus);
	  pParticle(PiMinus,production,piplus);
	  if (printGamma) {
	    pParticle(Gamma,production,gamma1);
	    pParticle(Gamma,production,gamma2);
	  }
	  else {
	    pParticle(Eta,production,eta);
	  }
	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}


/*------------------- test100 ---------------------------------------*/



void test100(int argc, char *argv[])
{
#define CONV 57.295791433

  int icount = 0;
  int printGamma = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double 
    masslow = 3 * PI_MASS,
    masshigh = 0.,
    t_max,
    slope = 10.0,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    virtualBeam,
    W,
    electron,
    resonance,
    recoil,
    piminus,
    piplus,
    pizero,
    isobar1;
  lorentzTransform Boost;

  fourVec CMbeam;

  
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass = BEAM_MASS;
  int printBaryon = 0;

  double ep_low = 1,ep_high = 4,theta_low = 0.5,theta_high = 1.2;
  double ep,thet,thetg,wsq,qx,qy;
  double px,py,pz,phi;

  float tMin;
  float setMassHigh = 0.0;

  int debug = 0;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof (ptr);
	break;
      case 'U':
	ptr++;
	setMassHigh = atof (ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'D':
	debug = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'e':
	ep_low = atof(++ptr);
	break;
      case 'E':
	ep_high = atof(++ptr);
	break;
      case 't':
	theta_low = atof(++ptr);
	break;
      case 'T':
	theta_high = atof(++ptr);
	break;

      case 'h':
	UsageM100(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double pi_pi_pi_threshold = 3 * PI_MASS;
  double pi_pi_threshold = 2 * PI_MASS;
  if(masslow < pi_pi_pi_threshold)
    masslow = pi_pi_pi_threshold;

 

  while (maxevents) {

    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);


    /*
     *-- beam and target in lab frame
     */

    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow (ELEC_MASS, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /* get an ep and a wsq */

    do {
      ep = randm(ep_low,ep_high);
      thetg = randm(theta_low,theta_high);
      thet = thetg/CONV;   // convert to radians
      qx = Qsq(beam.t(),ep,thet);
      wsq = Wsq(beam.t(),ep,thet);
    } while (1 < 0);
    phi = randm(0,2.0 * M_PI);
    pz = cos(thet) * ep;
    px = sin(thet) * cos(phi) * ep;
    py = sin(thet) * sin(phi) * ep;
      
    electron = fourVec( sqrt(ep * ep + pow (ELEC_MASS, 2.0)),threeVec(px,py,pz));

    // the W system is left- W in lab

    W = beam + target  - electron;
    qy = -(beam - electron).lenSq();
 
    double CMenergy = sqrt(wsq);

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */

    masshigh = (setMassHigh < pi_pi_pi_threshold) ? CMenergy - PROTON_MASS : setMassHigh;

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < pi_pi_pi_threshold);
    isobar1_mass = randm(pi_pi_threshold, (resonance_mass - PI_MASS));

    double resonance_p = CMmomentum (CMenergy, resonance_mass, PROTON_MASS);
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    // Only isotropic right now
    costheta = randm(-1.0,1.0);


    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (PROTON_MASS, 2.0)),zeroVec - resonance.V());

    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, PI_MASS);
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

    isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    piminus.set(sqrt(isobar1.V().lenSq() + pow (PI_MASS, 2.0)),zeroVec -isobar1.V());




    /*
     *  now do decay in isobar1 rest frame
     */
    // pi pi case
    double piplus_p = CMmomentum (isobar1_mass, PI_MASS,PI_MASS);

    piplus.polar( piplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   piplus.polar( piplus_p,acos (1.0),0.0);
    piplus.t(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)));
    pizero.set(sqrt (piplus.V().lenSq () + pow (PI_MASS, 2.0)),zeroVec - piplus.V());


    /*
     *  compute lorentz factor
     */
    LorentzFactor = (resonance_p * isobar1_p * piplus_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

        // boost from pi_pi(rho) rest frame to pi_pi_pi(a2/M) rest frame

        tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        piplus = Boost * piplus;
        pizero = Boost * pizero;


	// boost from pipipi to resonance

	tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);
	piplus = Boost * piplus;
	pizero = Boost * pizero;
	piminus = Boost * piminus;


        // boost to W 
        tmp.set(W.t(),zeroVec - W.V());
	Boost.set(tmp);

        isobar1 = Boost * isobar1;
        piplus = Boost * piplus;
        pizero = Boost * pizero;
        piminus = Boost * piminus;
	recoil = Boost * recoil;

	// boost from W to lab

	Boost.set(target);
	W *= Boost;
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        isobar1 = Boost * isobar1;
        piplus = Boost * piplus;
        pizero = Boost * pizero;
        piminus = Boost * piminus;

  
        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

	if (debug) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "W:\n  ";
	  W.print();
	  std::cerr << "Beam  in CM:\n  " ;
	  CMbeam.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "isobar1\n";
	  std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	  std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	  std::cerr << "  isobar1:\n    ";
	  isobar1.print();
	  std::cerr << "pi+ :\n  ";
	  piplus.print();
	  std::cerr << "pi-2:\n  ";
	  pizero.print();
	  std::cerr << "pi-1:\n  ";
	  piminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cerr << pow (~(pizero + piplus), 2.0) << " "
	       << pow (~(pizero + piminus), 2.0) << " "
	       << pow (~(piplus + piminus), 2.0) << std::endl;
	}
	else if (txt2part_style) {
	  if (printGamma) {
	    std::cout << "6" << std::endl;
	  } 
	  else { 
	    std::cout << "5" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;
	  pParticle_txt2part(Electron,production,electron);
	  pParticle_txt2part(Proton, production, recoil);
	  pParticle_txt2part(PiMinus, production, piminus);
	  pParticle_txt2part(PiPlus, production, piplus);
	  if (printGamma){
	    //      pParticle_txt2part(Gamma, production, gamma1);
	    //  pParticle_txt2part(Gamma, production, gamma2);
	  }
	  else{
	    pParticle_txt2part(Pi0,production,pizero);
	  }
	}
	else if (gamp) {
	   std::cout << (printBaryon ? "6" : "5") << std::endl;
	   pParticleGamp(Electron,beam);
	   pParticleGamp(Electron,electron);
	   pParticleGamp(PiPlus,piplus);
	   pParticleGamp(PiMinus,piminus);
	   pParticleGamp(Pi0,pizero);
	   if (printBaryon)
	     pParticleGamp(Proton,recoil);
	}
	else {
	  /*
	   *  write event
	   */
	  pParticle(Electron,production,electron);
	  pParticle(PiMinus,production,piplus);
	  pParticle(PiPlus,production,piminus);
	  pParticle(Pi0,production,pizero);
	  if (printBaryon)
	    pParticle(Proton,production,recoil);
	}
	if (Print) {
	  std::cout << "DBG " << ep << " " << thetg << " " << qx << " " << qy << " " << wsq << " ";
	  std::cout << std::endl;
	}
	  
	maxevents--;
      }
    }
  }
}

void doubleDalitz (char *ProcessName)
{
    std::cerr << ProcessName << " test bed\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass \n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#        : beam mass\n";
    std::cerr << "\t-j#        : beam id\n";
    std::cerr << "\t-B         : Print baryon\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}
void UsagedoubleDalitz (char *ProcessName)
{
}

void UsageM100 (char *ProcessName)
{
    std::cerr << ProcessName << " test bed\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass \n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#        : beam mass\n";
    std::cerr << "\t-j#        : beam id\n";
    std::cerr << "\t-B         : Print baryon\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}



/*----------------End of specific final states -----------------------*/

int comp_double (const void* a, const void* b)
{

    double* da = (double*) a;
    double* db = (double*) b;

    if(*da<*db){
	return -1;
    }
    else if(*da>*db){
	return 1;
    }
    else {
	return 0;
    }
}

double randm (double low, double high)
{
  return ((high - low) * drand48 () + low);
}

double CMmomentum (double cm_engy, double m1, double m2)
{
    double A = cm_engy*cm_engy;
    double B = (m1 + m2)*(m1 + m2);
    double C = (m1 - m2)*(m1-m2);
    double D = (A - B) * (A - C);
    if (D < 0)
      throw("CMmomentum illegal");
    return( ((D < 0 ) ? -sqrt(-D) : sqrt(D))/(2.0 * cm_engy));
}

void distribute_beam(vector4_t *beam)
{
 float r1,r2,r3,phi,theta,p,r;

 r1=gaussian_gen();
 beam->t=BEAM_EN_CENTER + r1*BEAM_EN_SIGMA;
 r2=gaussian_gen();
 phi=BEAM_PHI_CENTER + r2*BEAM_PHI_SIGMA;
 r3=gaussian_gen();
 theta=BEAM_THETA_CENTER + r3*BEAM_THETA_SIGMA;
 p=sqrt((beam->t*beam->t) - (PI_MASS*PI_MASS));
 r=p*sin(theta);
 beam->space.x=r*cos(phi);
 beam->space.y=r*sin(phi);
 beam->space.z=p*cos(theta);
}

void distribute_vertex(vector3_t *v)
{
 float r1,r2;

 r1=gaussian_gen();
 v->x=VERTEX_X_CENTER + r1*VERTEX_X_SIGMA;
 r2=gaussian_gen();
 v->y=VERTEX_Y_CENTER + r2*VERTEX_Y_SIGMA;
 v->z=randm(targetZ0,targetZ1);
}

float gaussian_gen()
{
 double r1,r2=9999,g1=1;

 while (r2>g1)
   {
    r1=randm(-3,3);
    g1=exp(-1*(r1*r1)/2)/sqrt(2.*M_PI);
    r2=randm(0,1);
   }
 return(r1);
}

void Usage (char *ProcessName)
{
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t\t: print dalitz plot masses\n";
    std::cerr << "\t-p\t\t: print events\n";
    std::cerr << "\t-ofile\t: write to file\n";
    std::cerr << "\t-h\t: print useage\n";
    std::cerr << flush;
}

void UsageM1 (char *ProcessName)
{
    std::cerr << ProcessName << " generates 3 pion monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e         : use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#        : beam mass\n";
    std::cerr << "\t-j#        : beam id\n";
    std::cerr << "\t-B         : Print baryon\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}

void UsageM10 (char *ProcessName)
{
    std::cerr << ProcessName << " generates 3 pion monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e         : use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#        : beam mass\n";
    std::cerr << "\t-B         : Print baryon\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}

void GeneralUsage()
{
  std::cerr << "\t-P#\t\tSet primary beam momentum to # (default: 4.0 GeV)\n";
  std::cerr << "\t-E#\t\tSet primary electron beam momentum to # (default: 4.0 GeV)\n";
  std::cerr << "\t-r#,#\t\tGenerate beam momentum uniform over given range" << std::endl;
  std::cerr << "\t-z[#,#]\t\tGenerate vertex z  over given range (default = -9.0,9.0)" << std::endl;
  std::cerr << "\t-b#,#\t\tGenerate a brem beam over given range" << std::endl;
  std::cerr << "\t-j#\t\tInput Beam id (default = Gamma [1])" << std::endl;
  std::cerr << "\t-I\t\tGenerate events isotropically" << std::endl;
  std::cerr << "\t-T\t\tMake events compatable with txt2part" << std::endl;
  std::cerr << "\t-G\t\tMake events compatable with gamp" << std::endl;
  std::cerr << "\t-d\t\tPrint dalitz plot variables" << std::endl;
  std::cerr << "\t-v\t\tverbose mode" << std::endl;
  std::cerr << "\t-A\t\tprint all particles" << std::endl;
  std::cerr << "\t-t#\t\texponential slope" << std::endl;
}

void MUsage(char *ProcessName)
{
  std::cerr << ProcessName << " generates phase space monte carlo events\n\n\n";
  std::cerr << "First you must choose a mode: -M# \n";
  std::cerr << "\t1:\t\tep -> e' pi+ pi- pi- p" << std::endl;
  std::cerr << "\t2:\t\tep -> e' pi+ pi- n" << std::endl;
  std::cerr << "\t3:\t\tep -> e' n pi+" << std::endl;  
  std::cerr << "\t4:\t\tep -> e' p pi0" << std::endl;
  std::cerr << "\t5:\t\tgamma p -> K+ Ks pi- p" << std::endl;
  std::cerr << "\t6:\t\tgamma p -> n pi+" << std::endl;
  std::cerr << "\t7:\t\tgamma p -> pi+ pi+ pi- n" << std::endl;
  std::cerr << "\t8:\t\tbeam p -> pi+ pi- pi0 p " << std::endl;
  std::cerr << "\t9:\t\tgamma p -> pi0 p " << std::endl;
  std::cerr << "\t10:\t\tep -> e' pi+ pi- pi0 p" << std::endl;
  std::cerr << "\t11:\t\tgamma  p ->(pi+ pi-) p" << std::endl;
  std::cerr << "\t12:\t\tgamma  p ->pi+/- (pi-/+ p)" << std::endl;
  std::cerr << "\t13:\t\te p -> e' pi+ pi- pi0 p" << std::endl;
  std::cerr << "\t14:\t\tpi p -> K+ K- n" << std::endl; 
  std::cerr << "\t15:\t\tgamma p -> e+ e- p" << std::endl;
  std::cerr << "\t16:\t\tbeam + p -> omega pi+ pi- (p or n) (default beam = gamma)" << std::endl;
  std::cerr << "\t17:\t\tbeam + p -> omega pi+ pi0 (p or n) (default beam = gamma)" << std::endl;
  std::cerr << "\t18:\t\tbeam + p -> K+ K- pi0 (p or n) (default beam = gamma)" << std::endl;
  std::cerr << "\t19:\t\tbeam + p -> K+ K- eta (p or n) (default beam = gamma)" << std::endl;
  std::cerr << "\t20:\t\tbeam + p -> K+ (K- p) (default beam = gamma)" << std::endl;
  std::cerr << "\t21:\t\tbeam p -> pi+ pi- eta p " << std::endl; 
  std::cerr << "\t22:\t\tbeam p -> pi+ pi- p " << std::endl;
  std::cerr << "\t23:\t\tbeam p -> pi+ pi+ pi-  (p,n) " << std::endl; 
  std::cerr << "\t24:\t\tbeam p -> K+ K- p " << std::endl;
  std::cerr << "\t25:\t\tbeam p -> (pbar p) p " << std::endl;
  std::cerr << "\t26:\t\tbeam p -> eta' pi p" << std::endl;
  std::cerr << "\t27:\t\tbeam p -> pi+ pi- pi-  (p,n) " << std::endl; 
  std::cerr << "\t28:\t\tbeam p -> pbar p p 3-body phase space" << std::endl;
  std::cerr << "\t29:\t\tbeam p -> K+ K- pi+  (p,n) " << std::endl; 
  std::cerr << "\t30:\t\tbeam p -> phi eta  (p,n) " << std::endl; 
  std::cerr << "\t31:\t\tbeam p -> K+ (K- p) " << std::endl;
  std::cerr << "\t32:\t\tbeam p -> pi+ (pi- p) " << std::endl;
  std::cerr << "\t33:\t\tbeam p -> pi- (pi+ p) " << std::endl; 
  std::cerr << "\t34:\t\tbeam p -> K+ Kl n " << std::endl;
  std::cerr << "\t35:\t\tbeam p -> p pi+ pi-  3-body phase space" << std::endl; 
  std::cerr << "\t36:\t\tbeam p -> Kshort (-> pi+ pi-) K+ n" << std::endl; 
  std::cerr << "\t37:\t\tbeam p -> e+ e- p (exp distribution of m(e+e-)" << std::endl;
  std::cerr << "\t38:\t\tbeam p -> pi0 pi0  p (Double Dalitz)" << std::endl;
  std::cerr << "\t39:\t\tbeam p -> p pi0 pi0  3-body phase space" << std::endl;  
  std::cerr << "\t40:\t\tbeam He4 -> alpha pi0 eta " << std::endl;   
  std::cerr << "\t41:\t\tbeam p -> p pi0 eta " << std::endl;   
   std::cerr << "\t42:\t\tbeam p -> Kshort (-> pi+ pi-) KL p" << std::endl;
    std::cerr << "\t43:\t\tbeam p -> phi omega p" << std::endl;
  std::cerr << "\nAfter choosing the appropriate mode, call ppgen again with the -h option, eg:";
  std::cerr << "\nppgen -M# -h\nThis will document options appropriate to that particular mode\n"; 
}


void UsageM2 (char *ProcessName)
{
    std::cerr << ProcessName << " generates 2 pion monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e         : use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}

void UsageM3 (char *ProcessName)
{
    std::cerr << ProcessName << " generates neutron pi+  monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e         : use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}



void UsageM5 (char *ProcessName)
{
    std::cerr << ProcessName << " generates k+ ks pi- monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e         : use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}
void UsageM6(char *ProcessName)
{
    std::cerr << ProcessName << " generates gamma p -> neutron pi+  monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-e         : use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}


void UsageM7 (char *ProcessName)
{
    std::cerr << ProcessName << " generates 3 pion monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e         : use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#        : beam mass\n";
    std::cerr << "\t-B         : Print baryon\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}
void UsageM8 (char *ProcessName)
{
    std::cerr << ProcessName << " generates 3 pion monte carlo events (pi+ pi- pi0)\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#\t: beam mass\n";
    std::cerr << "\t-B\t: Print baryon\n";
    std::cerr << "\t-g\t: Print gammas (default = print pi0)" << std::endl;
    std::cerr << "\t-p\t: print events\n";
    std::cerr << "\t-h\t: print usage\n";
    std::cerr << "\t\tdefault slope: 10.0" << std::endl;
    std::cerr << flush;
    
}
void UsageM9 (char *ProcessName)
{
    std::cerr << ProcessName << " generates p pi0  pion monte carlo events \n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#\t: beam mass\n";
    std::cerr << "\t-B\t: Print baryon\n";
    std::cerr << "\t-g\t: Print gammas (default = print pi0)" << std::endl;
    std::cerr << "\t-p\t: print events\n";
    std::cerr << "\t-h\t: print useage\n";
    std::cerr << flush;
    
}
void UsageM11(char *ProcessName)
{
    std::cerr << ProcessName << " generates 2 pion monte carlo events- photon beam (pi+ pi-)  - pi+ pi- system forward\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#\t: beam mass\n";
    std::cerr << "\t-B\t: Print baryon\n";
     std::cerr << "\t-p\t: print events\n";
    std::cerr << "\t-h\t: print useage\n";
    std::cerr << "\t\tdefault slope: 10.0" << std::endl;
    std::cerr << flush;
    
}
void UsageM12(char *ProcessName)
{
    std::cerr << ProcessName << " generates 2 pion monte carlo events- photon beam (pi+ pi-)  - pi  forward\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#\t: beam mass\n";
    std::cerr << "\t-B\t: Print baryon\n";
     std::cerr << "\t-p\t: print events\n";
     std::cerr << "\t=cX\t: charge of the forward pi, (+/-) (default = -)" << std::endl;
    std::cerr << "\t-h\t: print useage\n";
    std::cerr << "\t\tdefault slope: 10.0" << std::endl;
    std::cerr << flush;
    
}

void UsageM13(char *ProcessName)
{
  std::cerr << ProcessName << " e p -> e' p pi+ pi- pi0\nDEVELOPMENT- DO NOT USE!!!!\n\n";
  std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
  GeneralUsage();
  std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
  std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
  std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
  std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
  std::cerr << "\t-e\t: use 1/E factor in LIPS factor\n";
  std::cerr << "\t-d\t: print dalitz plot masses (to stdout)\n";
  std::cerr << "\t-b#\t: beam mass\n";
  std::cerr << "\t-B\t: Print baryon\n";
  std::cerr << "\t-q#\t:Q^2 lower limit" << std::endl;
  std::cerr << "\t-Q#\t:Q^2 upper limit" << std::endl;
  std::cerr << "\t-w#\t:W lower limit" << std::endl;
  std::cerr << "\t-W#\t:W upper limit" << std::endl;
  std::cerr << "\t-p\t: print events\n";
  std::cerr << "\t-h\t: print useage\n";
  std::cerr << "\t\tdefault slope: 10.0" << std::endl;
  std::cerr << flush;
    
}

void UsageM14 (char *ProcessName)
{
    std::cerr << ProcessName << " generates b p -> K+ k- n monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents : use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass  : lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass  : upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e         : use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d         : print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-p         : print events\n";
    std::cerr << "\t-ofile     : write to file\n";
    std::cerr << "\t-h         : print useage\n";
    std::cerr << flush;
    
}

void UsageM15(char *ProcessName)
{
    std::cerr << ProcessName << " generates 1 e+ and 1 e- monte carlo events- photon beam (e+ e-)  - e+ e- system forward\n\n\n";
    std::cerr << "Usage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#\t\t: beam mass\n";
    std::cerr << "\t-B\t\t: Print baryon\n";
     std::cerr << "\t-p\t\t: print events\n";
    std::cerr << "\t-h\t\t: print useage\n";
    std::cerr << "\t\tdefault slope: 10.0" << std::endl;
    std::cerr << flush;
    
}
void UsageM16(char *ProcessName)
{
    std::cerr << ProcessName << " generates omega pi pi final state\n\n\n";
    std::cerr << "Usage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-B\t\t: Print baryon\n";
     std::cerr << "\t-p\t\t: print events\n";
    std::cerr << "\t-h\t\t: print useage\n";
    std::cerr << "\t\tdefault slope: 4.0" << std::endl;
    std::cerr << flush;
    
}
void UsageM26(char *ProcessName)
{
    std::cerr << ProcessName << " generates etaprime pi final state\n\n\n";
    std::cerr << "Usage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-B\t\t: Print baryon\n";
     std::cerr << "\t-p\t\t: print events\n";
    std::cerr << "\t-h\t\t: print useage\n";
    std::cerr << "\t\tdefault slope: 4.0" << std::endl;
    std::cerr << flush;
    
}

void Usagecpcmn0N(char *ProcessName)
{
    std::cerr << ProcessName << " generates beam + proton -> C+ C- N0  monte carlo events (ex: K+ K- pi0)\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-b#\t: beam mass\n";
    std::cerr << "\t-B\t: Print baryon\n";
    std::cerr << "\t-g\t: Print gammas (default = print pi0)" << std::endl;
    std::cerr << "\t-p\t: print events\n";
    std::cerr << "\t-h\t: print useage\n";
    std::cerr << "\t\tdefault slope: 10.0" << std::endl;
    std::cerr << flush;
    
}
void Usagecpcm(char *ProcessName)
{
    std::cerr << ProcessName << " generates beam + proton -> C+ C- recoil  monte carlo events (ex: K+ K-)\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-B\t: Print baryon\n";
    std::cerr << "\t-p\t: print events\n";
    std::cerr << "\t-h\t: print useage\n";
    std::cerr << "\t\tdefault slope: 10.0" << std::endl;
    std::cerr << flush;
    
}
void Usagec1c2(char *ProcessName)
{
    std::cerr << ProcessName << " generates beam + proton -> C1 C2  (C1 -> a b, C2 -> c,d) recoil  monte carlo events\n\n\n";
    std::cerr << "Useage: " << ProcessName << " -mmaxevents -llfevents -Lmass -Umass -D -p -ofile -h\n";
    GeneralUsage();
    std::cerr << "\t-mmaxevents\t: write maxevents, default 999999\n";
    std::cerr << "\t-llfevents\t: use lfevents to find max. lorentz factor, default 3000\n";
    std::cerr << "\t-Llowmass\t: lower bound for X mass, default 1.5 GeV\n";
    std::cerr << "\t-Ulowmass\t: upper bound for X mass, default 2.5 GeV\n";
    std::cerr << "\t-e\t: use 1/E factor in LIPS factor\n";
    std::cerr << "\t-d\t: print dalitz plot masses (to stdout)\n";
    std::cerr << "\t-B\t: Print baryon\n";
    std::cerr << "\t-p\t: print events\n";
    std::cerr << "\t-h\t: print useage\n";
    std::cerr << "\t\tdefault slope: 10.0" << std::endl;
    std::cerr << flush;
    
}



void generateBeamMCinsideTarget(vector3_t *v,vector3_t *p)
{
  if (UseTarget) {
    distribute_vertex(v);
  }
  else {
    v->x = v->y = 0.0;
    v->z = z;
  }

  p->x = p->y = 0.0; 
  p->z = GetBeamP();
}


void pParticle_txt2part(Particle_t type,threeVec v,fourVec p)
{
  std::cout << (int)type << std::endl;;
  std::cout << p.t() << " " << p.x() << " " << p.y() << " " << p.z() << std::endl;
  std::cout << v.x() << " " << v.y() << " " << v.z() << std::endl;; 
}
void pParticleGamp(Particle_t type,fourVec p)
{
  std::cout << (int)type << " ";
  std::cout << (int) Q(type) << " ";
  std::cout << p.x() << " " << p.y() << " " << p.z() << " " << p.t() << std::endl;
}
void pParticle(Particle_t type,threeVec v,fourVec p)
{

  std::cout << (int) type << " ";
  std::cout <<   v.x() << " " << v.y() << " " << v.z() << " ";
  std::cout << p.t() << " " << p.x() << " " << p.y() << " " << p.z() << std::endl;

}

double GetBeamP()
{
  double ret;
    if (UseRange) {
    ret = randm(range0,range1);
    }
  else if (UseBrem) {
    ret = brem(brem0,brem1); 
  }
  else
    ret = pbeamZ;

 
  return(ret);
  
}

double brem(double brem0,double brem1)
{
  double r = randm(0.,1.);
  double ret = brem0 * exp(r *  log(brem1/brem0));
  
  return(ret);

}

float e(float p,float m)
{
  return(sqrt(p * p + m * m));
}
float s(float plab,float m1,float m2)
{

  return(m1*m1 + m2*m2 + 2.0 * sqrt(plab * plab + m1 * m1) * m2);
}
float wcm(float p,float m1, float m2)
{
  return(sqrt(p*p + m1*m1) + sqrt(p*p + m2 * m2));
}
float pprime(float w, float m, float m3)
{
  return( sqrt( pow( (w*w - m*m + m3*m3)/(2.0 * w),2) - m3*m3));
}

double getT(double tMin,double slope)
{
  double t;
  int icount = 0;
  do {
    if((icount++) > 100000) {
      std::cerr << "getT: exiting... too many counts\t tMin: " << tMin << "\tslope: " << slope  << std::endl;
      exit(0);
    }

    t = log(randm(0.,1.0))/slope;
  } while( t > tMin);
  return(t);
}
double tprimeExp(double plab,double beamMass,double targetMass ,double resonanceMass,double recoilMass,double slope) {
        float tMin = tmin(plab,beamMass,targetMass,resonanceMass,recoilMass);
        float tMax = tmax(plab, beamMass, targetMass, resonanceMass, recoilMass);
        float W = sqrt(s(plab, beamMass, targetMass));
        float resonance_p = pprime(W, resonanceMass, recoilMass);
        float beam_p = pprime(W,beamMass,targetMass);
        double t0 = fabs(tMin);
        double t1 = fabs(tMax);
        double A = slope/(exp(-slope * t0) - exp(-slope * t1));
        double r = randm(0.0,1.0);
        double arg = -slope * r/A + exp(-slope * t0);
        double ret = -log(arg)/slope;
        return(-ret);
}
double cosThetaExp(double plab,double beamMass,double targetMass ,double resonanceMass,double recoilMass,double slope) {
        double tt = tprimeExp(plab,beamMass,targetMass,resonanceMass,recoilMass,slope);
        double tMin = tmin(plab,beamMass,targetMass,resonanceMass,recoilMass);
        double tMax = tmax(plab,beamMass,targetMass,resonanceMass,recoilMass);
        double tx = tMin - tt;
        float W = sqrt(s(plab, beamMass, targetMass));
        float resonance_p = pprime(W, resonanceMass, recoilMass);
        float beam_p = pprime(W,beamMass,targetMass);
        double ret = 1.0 - (tx/(2.0 * resonance_p * beam_p));
        if (fabs(ret) > 1.0) {
            cerr << "ERROR " << tt << " " << tx << " " << tMin <<  " " << tMax << " " << W << " " << resonance_p << " " << beam_p << " " << ret << " " << resonanceMass << " " << slope << endl;
            ret = 1.0;
        }
        return ret;
    }


float tmin(float plab,float ma,float mb,float mc,float md)
{
  float w = sqrt(s(plab,ma,mb));
  float pcma = pprime(w,ma,mb);
  float pcmc = pprime(w,mc,md);
  float ea = e(pcma,ma);
  float ec = e(pcmc,mc);
  return(ma*ma + mc*mc - 2.0 * ea * ec + 2.0 * pcma * pcmc);
}
float tmax(float plab,float ma,float mb,float mc,float md)
{
  float w = sqrt(s(plab,ma,mb));
  float pcma = pprime(w,ma,mb);
  float pcmc = pprime(w,mc,md);
  float ea = e(pcma,ma);
  float ec = e(pcmc,mc);
  return(ma*ma + mc*mc - 2.0 * ea * ec - 2.0 * pcma * pcmc);
}


float theta(float QSQ,float W,float E) 
{

  double k;

  k = (2. * PROTON_MASS * QSQ)/((2.0 * PROTON_MASS * E + PROTON_MASS * PROTON_MASS - W * W - QSQ) * (4.0 * E));
  if (fabs(k) <= 1.0 && k > 0.0) {
    return(2.0 * asin(sqrt(k)));
  }
  else
    return(-1000.0);
}

float eprime(float theta,float Qsq,float E)
{

  return(Qsq/(4. *  E * sin(theta/2.0) * sin(theta/2.0)));
}


double Qsq(double E,double Ep,double theta)
{
  return(4.0 * E * Ep * sin(theta/2.0) * sin(theta/2.0));
}

double nQsq(double E,double Ep,double theta)
{
  double p = sqrt(E*E - ELECTRON_MASS * ELECTRON_MASS);
  double pp = sqrt(Ep*Ep - ELECTRON_MASS * ELECTRON_MASS);
  double ret = 2.0 * (ELECTRON_MASS * ELECTRON_MASS - E * Ep + p * pp * cos(theta));
  return(ret);
}
	 

double Wsq(double E,double Ep,double theta)
{
  return(2.0 * PROTON_MASS * (E - Ep) - 2.0 * E * Ep * (1.0 - cos(theta)) + PROTON_MASS * PROTON_MASS);
}


  
int  Q(Particle_t  pid)
{
  int ret = 0;
  switch (pid) {
  case Alpha:
    ret = 2;
    break;
  case PiPlus:
  case Proton:
  case KPlus:
  case Positron:
    ret = 1;
    break;
  case PiMinus:
  case Electron:
  case KMinus:
  case AntiProton:
    ret = -1;
    break;
  case Gamma:
  case Neutron:
  case AntiNeutron:
  case KLong:
  case KShort:
  case Pi0:
  case Rho0:
  case Eta:
  case EtaPrime:
  case omega:
  case phiMeson:
    ret = 0;
    break;
  default:
    ret =-100;
    break;
  }
  return(ret);
}
double Mass(Particle_t  pid)
{
#define ALPHA_MASS 3.73
  double  ret = 0;
  switch (pid) {
  case Alpha:
    ret = ALPHA_MASS;
    break;
  case phiMeson:
    ret = PHI_MASS;
    break;
  case Eta:
    ret = ETA_MASS;
    break;
  case PiPlus:
  case PiMinus:
    ret = PI_MASS;
    break;
  case Proton:
  case AntiProton:
    ret = PROTON_MASS;
    break;
  case Neutron:
  case AntiNeutron:
    ret = NEUTRON_MASS;
    break;

  case KPlus:
  case KMinus:
    ret = KCHARGED_MASS;
    break;

  case Positron:
  case Electron:
    ret = 0.00051;
    break;
  case Gamma:
    ret = 0.0;
    break;
  case KLong:
  case KShort:
    ret = KZERO_MASS;
    break;
  case Pi0:
    ret = PI0_MASS;
    break;
  case omega:
    ret = OMEGA_MASS;
    break;
    
  default:
    ret =0.0;
    break;
  }
  return(ret);
}

void PrintParticleID()
{
        std::cerr << " Unknown = 0" << std::endl;
        std::cerr << " Gamma = 1" << std::endl;
        std::cerr << " Positron = 2" << std::endl;
        std::cerr << " Electron = 3" << std::endl;
        std::cerr << " Neutrino = 4" << std::endl;
        std::cerr << " MuonPlus = 5" << std::endl;
        std::cerr << " MuonMinus = 6" << std::endl;
        std::cerr << " Pi0 = 7" << std::endl;
        std::cerr << " PiPlus = 8" << std::endl;
        std::cerr << " PiMinus = 9" << std::endl;
        std::cerr << " KLong = 10" << std::endl;
        std::cerr << " KPlus = 11" << std::endl;
        std::cerr << " KMinus = 12" << std::endl;
        std::cerr << " Neutron = 13" << std::endl;
        std::cerr << " Proton = 14" << std::endl;
        std::cerr << " AntiProton = 15" << std::endl;
        std::cerr << " KShort = 16" << std::endl;
        std::cerr << " Eta = 17" << std::endl;
        std::cerr << " Lambda = 18" << std::endl;
        std::cerr << " SigmaPlus = 19" << std::endl;
        std::cerr << " Sigma0 = 20" << std::endl;
        std::cerr << " SigmaMinus = 21" << std::endl;
        std::cerr << " Xi0 = 22" << std::endl;
        std::cerr << " XiMinus = 23" << std::endl;
        std::cerr << " OmegaMinus = 24" << std::endl;
        std::cerr << " AntiNeutron = 25" << std::endl;
        std::cerr << " AntiLambda = 26" << std::endl;
        std::cerr << " AntiSigmaMinus = 27" << std::endl;
        std::cerr << " AntiSigma0 = 28" << std::endl;
        std::cerr << " AntiSigmaPlus = 29" << std::endl;
        std::cerr << " AntiXi0 = 30" << std::endl;
        std::cerr << " AntiXiPlus = 31" << std::endl;
        std::cerr << " AntiOmegaPlus = 32" << std::endl;
        std::cerr << " Deuteron = 45" << std::endl;
        std::cerr << " Triton = 49" << std::endl;
        std::cerr << " Rho0 = 57" << std::endl;
        std::cerr << " RhoPlus = 58" << std::endl;
        std::cerr << " RhoMinus = 59" << std::endl;
        std::cerr << " omega = 60" << std::endl;
        std::cerr << " EtaPrime = 61" << std::endl;
        std::cerr << " phiMeson = 62" << std::endl;
}





//--------------------- New code -------------------

// beam p -> c1 c2 (p,n)  c1-> p1_1 p1_2  c2 -> p2_1 p2_2

void c1c2(int argc,char **argv,Particle_t Beam,Particle_t Target,Particle_t C1,Particle_t P1_1,Particle_t P1_2,Particle_t C2,Particle_t P2_1,Particle_t P2_2,Particle_t Recoil)
{
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  Particle_t Baryon;
  double 
    masslow = Mass(C1) + Mass(C2),
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    c1,
    c2,
    p1_1,p1_2,p2_1,p2_2;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass;
  float baryonmass;
  int printBaryon = 0;
  int debug = 0;  
  double targetMass;
  int massHighSet = 0;


  /* generate vertices */
  threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
        massHighSet = 1;
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'D':
	debug = 1;
	break;
      case 'h':
	Usagec1c2(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 10.0;
  int qtot = Q(Beam) + Q(Proton) - Q(C1) - Q(C2);
  double c1c2Threshold =  Mass(C1) + Mass(C2);

  Baryon = Recoil;

  beamMass = Mass(Beam);
  baryonmass = Mass(Baryon);
  targetMass = Mass(Target);

  if(masslow < c1c2Threshold)
    masslow = c1c2Threshold;

 

  while (maxevents) {
      float thisMassHigh = masshigh;
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (targetMass,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(!massHighSet)
      thisMassHigh = CMenergy - Mass(Baryon);
    
    float availableMass = CMenergy - Mass(Baryon);
    
    if (thisMassHigh > availableMass) {
        thisMassHigh = availableMass;
    }
    
    
    if (masshigh < c1c2Threshold) {
      std::cerr << "Meson high mass below 2 particle threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, thisMassHigh);
    } while (resonance_mass < c1c2Threshold);

    double beam_p      = CMmomentum (CMenergy, Mass(Beam), targetMass);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, Mass(Baryon));
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double t;

    double lf;

    tMin = tmin(pbeam.z,Mass(Beam),targetMass,resonance_mass,Mass(Baryon));

    if (Isotropic)
      costheta = randm(-1.0,1.0);
    else {

      t = getT(tMin,slope);
      costheta = 1 + (t - tMin)/(2 * resonance_p * beam_p);

    }


    resonance.polar(resonance_p,acos (costheta), randm (-M_PI,M_PI));
    resonance.t(resonance_E);


    /*
     *-- recoil particle
     */
    recoil.set(sqrt(resonance.V().lenSq() +  pow (Mass(Baryon), 2.0)),zeroVec - resonance.V());


    /*
     *  now do decay in isobar1 rest frame
     */
    

    lf = Decay2(resonance,~resonance,c1,Mass(C1),c2,Mass(C2));

    // c1 and c2 now in overall CM  rest frame

    lf *= Decay2(c1,Mass(C1),p1_1,Mass(P1_1),p1_2,Mass(P1_2));
    lf *= Decay2(c2,Mass(C2),p2_1,Mass(P2_1),p2_2,Mass(P2_2));

    // Boost everyone to lab frame
    Boost.set (CMtarget);
        resonance = Boost * resonance;
        recoil = Boost * recoil;
        c1 = Boost * c1;
        c2 = Boost * c2;
	p1_1 *= Boost;
	p1_2 *= Boost;
	p2_1 *= Boost;
	p2_2 *= Boost;	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "c1 :\n  ";
	  c1.print();
	  std::cerr << "c2:\n  ";
	  c2.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factors: " << LorentzFactor << " " << lf << "\n";
	}
	// calculate masses for dalitz plots
        
        if (!std::isnan(recoil.t())) {
        
	if (dalitz) {
	  std::cout << resonance_mass << " ";
	  std::cout << ~(recoil + c1 + c2) << " ";
	  std::cout << pow (~(recoil + c1), 2.0) << " "
	       << pow (~(recoil + c2), 2.0) << " "
	       << pow (~(c1 + c2), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << std::endl;
	} 
	else if (debug) {
	  std::cout << "T " << (beam - c1 - c2).lenSq() << " " << (target - recoil).lenSq() << " " << t << " " << tMin << std::endl;
	}
	if (txt2part_style){
	  std::cout << "4" << std::endl;
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Recoil, production, recoil);
	  pParticle_txt2part(C1, production, c1);
	  pParticle_txt2part(C2, production, c2);
	}
	else if (gamp) {
	  event evt;
	  
	  particle gBeam(PDGtable.get(pid2name(Beam)),Q(Beam));
	  particle gTarget(PDGtable.get(pid2name(Target)),Q(Target));
	  particle gRecoil(PDGtable.get(pid2name(Recoil)),Q(Recoil));
	  particle gC1(PDGtable.get(pid2name(C1)),Q(C1));
	  particle gC2(PDGtable.get(pid2name(C2)),Q(C2));
	  particle gP1_1(PDGtable.get(pid2name(P1_1)),Q(P1_1));
	  particle gP1_2(PDGtable.get(pid2name(P1_2)),Q(P1_2));
	  particle gP2_1(PDGtable.get(pid2name(P2_1)),Q(P2_1));
	  particle gP2_2(PDGtable.get(pid2name(P2_2)),Q(P2_2));

	  gBeam.set4P(beam);
	  gTarget.set4P(target);
	  gRecoil.set4P(recoil);
	  gC1.set4P(c1);
	  gC2.set4P(c2);
	  gP1_1.set4P(p1_1);
	  gP1_2.set4P(p1_2);
	  gP2_1.set4P(p2_1);
	  gP2_2.set4P(p2_2);

	  evt.beam(gBeam);
	  evt.target(gTarget);
	  evt.addfinal(gRecoil);
	  if (printAll) {
	    evt.addfinal(gC1);
	    evt.addfinal(gC2);
	  }
	  evt.addfinal(gP1_1);
	  evt.addfinal(gP1_2);
	  evt.addfinal(gP2_1);
	  evt.addfinal(gP2_2);
	  
	  std::cout << evt;
	  
	  
	}
        }
	maxevents--;
  }
} 

/* ---------------- beam p -> c1 c2, c1 -> d1 d2 -------------------------------------*/

/*  example pi- p -> pi0 n*/

void twobody (int argc, char *argv[],Particle_t Beam,Particle_t Target,Particle_t C1,Particle_t C2,int decay,Particle_t D1,Particle_t D2)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  int printDecay = 0;
  double 
    masslow = Mass(C1) + Mass(C2),
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0;
  fourVec 
    beam,
    target,
    c1,
    c2,
    d1,
    d2;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass;
  int debug = 0;
  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'g':
	printDecay = 1;
	break; 
      case 'D':
	debug = 1;
	break;
      case 'h':
	Usagecpcmn0N(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 10.0;
  double Threshold = Mass(C1) + Mass(C2);
  double targetmass = Mass(Target);
  int qtot = Q(Beam) + Q(Proton) - Q(C1) - Q(C2);
  if (qtot) {
    std::cerr << "illegal charge combination:\ttotal charge  = " << qtot << std::endl;
    exit(0);
  }


  beamMass = Mass(Beam);

  if(masslow < Threshold)
    masslow = Threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (targetmass,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy;
    if (masshigh < Threshold) {
      std::cerr << "High mass below 2particle threshold" << std::endl;
      exit(1);
    }

    double beam_p      = CMmomentum (CMenergy, Mass(Beam), targetmass);
    double c1_p        = CMmomentum(CMenergy,Mass(C1),Mass(C2));
    double c1_E        = sqrt (pow (Mass(C1), 2.0) + pow (c1_p, 2.0));
    double costheta;
    double t;

    tMin = tmin(pbeam.z,Mass(Beam),targetmass,Mass(C1),Mass(C2));

    t_max = 4. * beam_p * c1_p;
    expt_max = exp(-1.)/slope;


    if (Isotropic)
      costheta = randm(-1.0,1.0);

    else {
      t = getT(tMin,slope);
      costheta = 1 + (t - tMin)/(2 * c1_p * beam_p);
    }

    c1.polar(c1_p,acos (costheta), randm (-M_PI,M_PI));
    c1.t(c1_E);


    // recoil particle

    c2.set(sqrt(c1.V().lenSq() +  pow (Mass(C2), 2.0)),zeroVec - c1.V());

    
    // decay c1, in its rest frame
    if (decay) {
      double d1_p = CMmomentum (Mass(C1), Mass(D1),Mass(D2));
      d1.polar(d1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
      d1.t(d1_p);
      d2.set(d1_p,threeVec (0.0, 0.0, 0.0) - d1.V());
      d2.t(d1_p);
    }
    /*
     *  compute lorentz factor
     */
    LorentzFactor = FlatMass ? 1 : (c1_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      /* generate vertices */
      threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

	// boost d1,d2 to c1rest frame
	if (decay) {
	  tmp.set(c1.t(),zeroVec - c1.V());
	  Boost.set(tmp);
	  d1 = Boost * d1;
	  d2 = Boost * d2;
	}

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        c1 = Boost * c1;
        c2 = Boost * c2;
	if (decay) {
	  d1 = Boost * d1;
	  d2 = Boost * d2;
	}

	if (Print) {

	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "c1 :\n  ";
	  c1.print();
	  std::cerr << "c2:\n  ";
	  c2.print();
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  ;
	} 
	else if (debug) {
	  std::cout << "Z " << production.x() << " " << production.y() << " " << production.z() << std::endl;
	}
	else if (txt2part_style){
	  int ntxt;
	  if (printAll) {
	    ntxt = 4;
	  }
	  else if (decay && printDecay)
	    ntxt = 3;
	  else
	    ntxt = 2;

	  std::cout << ntxt << std::endl;
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(C2, production, c2);

	  if (decay && printDecay){
	    pParticle_txt2part(D1, production, d1);
	    pParticle_txt2part(D2, production, d2);
	  }
	  else{
	    pParticle_txt2part(C1, production, c1);
	  }
	}
	else if (gamp) {

	  event evt;
	  
	  particle gBeam(PDGtable.get(pid2name(Beam)),Q(Beam));
	  particle gTarget(PDGtable.get(pid2name(Target)),Q(Target));
	  particle gC2(PDGtable.get(pid2name(C2)),Q(C2));
	  particle gC1(PDGtable.get(pid2name(C1)),Q(C1));
	  particle gD1(PDGtable.get(pid2name(D2)),Q(D1));
	  particle gD2(PDGtable.get(pid2name(D2)),Q(D1));

	  gBeam.set4P(beam);
	  gTarget.set4P(target);
	  gC2.set4P(c1);
	  gC1.set4P(c2);
	  if (decay) {
	    gD1.set4P(d1);
	    gD2.set4P(d2);
	  }


	  evt.beam(gBeam);
	  evt.target(gTarget);
	  evt.addfinal(gC2);
	  if (decay && printDecay) {
	    evt.addfinal(gD1);
	    evt.addfinal(gD2);
	  }
	  else if (printAll) {
	    if (decay) {
	      evt.addfinal(gD1);
	      evt.addfinal(gD2);
	    }
	    evt.addfinal(gC1);
	  }
	  else {
	    evt.addfinal(gC1);
	  }
	  
	  std::cout << evt;
	}
	else {	  	  
	  /*
	   *  write event
	   */
	

	  if (printBeam)
	    pParticle(Beam,production,beam);

	  pParticle(C2,production,c2);

	  if (!decay || printAll)
	    pParticle(C1,production,c1);

	  if (decay && printDecay) {
	    pParticle(D1,production,d1);
	    pParticle(D2,production,d2);
	  }
	  else {
	    pParticle(C1,production,c1);
	  }
	}
	nevents++;
	if (verbose) {
	  if (!(nevents % 100)) 
	    std::cerr << nevents << "\r" << flush;
	} 
	maxevents--;
      }

    
    }
  }

  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

/*----------------End of twobody -----------------------------------*/

string pid2name(Particle_t  type) 
{
    switch (type) {
    case Alpha:
      return "alpha";
      break;
    case PiMinus:
    case PiPlus:
	return "pi";
	break;
    case KMinus:
    case KPlus:
	return "K";
	break;
    case Pi0:
	return "pi0";
	break;
    case Eta:
	return "eta";
	break;
    case phiMeson:
      return "phi(1020)";
      break;
    case omega:
	return ("omega(782)");
	break;
    case Proton:
	return "p";
	break;
	case AntiProton:
	return("pbar");
	break;
    case Neutron:
	return "n";
	break;
    case Gamma:
	return ("gamma");
	break;
    case Electron:
    case Positron:
	return ("e");
	break;
    case KShort:
    case KLong:
      return("K0");
      break;
    default:
	return "unknown";
	break;
    }
}

// --------------------  Decay2 -----------------------------------------------------

double Decay2(fourVec &Rest,double massRest,fourVec &P1,double massP1,fourVec &P2,double massP2)
{
  // decays Rest -> P1 P2 and returns phase space factor
  // fourVec P1 and P2 are in rest frame of Rest

  threeVec zeroVec(0.0,0.0,0.0);
  fourVec tmp;
  double p = CMmomentum(massRest,massP1,massP2);
  double costheta = randm(-1.0,1.0);
  double sintheta = sqrt(1 - costheta * costheta);
  double phi = randm(-M_PI,M_PI);
  lorentzTransform L(Rest);
  threeVec p1(p * cos(phi) * sintheta,p * sin(phi) * sintheta,p * costheta);
  threeVec p2 = zeroVec - p1;
  P1.set(sqrt(p * p + massP1 * massP1),p1);
  P2.set(sqrt(p * p + massP2 * massP2),p2);
  // Transform back

  tmp.set(Rest.t(),zeroVec - Rest.V());
  L.set(tmp);
  P1 *= L;
  P2 *= L;
  return(p);
}


  
  
  


// --------------------- 3 body phase space ------------------------------------------

void body3 (int argc, char *argv[],Particle_t Beam,Particle_t Cplus,Particle_t Cminus,Particle_t N0,int decay)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  int PrintGammas = 0;
  Particle_t Baryon;
  double 
    masslow = Mass(Cplus) + Mass(Cminus) + Mass(N0),
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    recoil,
    cminus,
    cplus,
    n0,
    gamma1,
    gamma2,
    isobar1;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass;
  float baryonmass;
  int printBaryon = 0;
  int printGamma = 0;
  int debug = 0;
  Particle_t Target = Proton;
  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'g':
	printGamma = 1;
	break; 
      case 'D':
	debug = 1;
	break;
      case 'h':
	Usagecpcmn0N(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 10.0;
  double cpcmn0Threshold = Mass(Cplus) + Mass(Cminus) + Mass(N0);
  double cpcmThreshold = Mass(Cplus) + Mass(Cminus); 
  int qtot = Q(Beam) + Q(Proton) - Q(Cplus) - Q(Cminus) -Q(N0);
  switch (qtot) {
  case 0:
    Baryon = Neutron;
    break;
  case 1:
    Baryon = Proton;
    break;
  default:
    std::cerr << "illegal charge combination:\ttotal charge  = " << qtot << std::endl;
    exit(0);
  }


  beamMass = Mass(Beam);

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);


    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();

    resonance_mass = ~(beam + target);
    isobar1_mass = randm(cpcmThreshold,resonance_mass - Mass(N0));

    double beam_p      = CMmomentum (CMenergy, Mass(Beam), PROTON_MASS);
    double resonance_p = 0.0;
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));



    resonance = fourVec(resonance_mass,threeVec(0.0,0.0,0.0));

    /*
     *  now do decay in resonance rest frame
     */
    double isobar1_p = CMmomentum (resonance_mass, isobar1_mass, Mass(N0));
     
    isobar1.polar( isobar1_p, acos (randm (-0.999999, 0.999999)), randm (-M_PI, M_PI)); 
    //   isobar1.polar( isobar1_p, acos (1.0), 0.0);

    isobar1.t(sqrt (isobar1.V().lenSq () + pow (isobar1_mass, 2.0)));
    
    n0.set(sqrt(isobar1.V().lenSq() + pow (Mass(N0), 2.0)),zeroVec -isobar1.V());




    /*
     *  now do decay in isobar1 rest frame
     */
    // c+ c-
    double cplus_p = CMmomentum (isobar1_mass, Mass(Cplus),Mass(Cminus));

    cplus.polar( cplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    //   cplus.polar( cplus_p,acos (1.0),0.0);
    cplus.t(sqrt (cplus.V().lenSq () + pow (Mass(Cplus), 2.0)));
    cminus.set(sqrt (cplus.V().lenSq () + pow (Mass(Cminus), 2.0)),zeroVec - cplus.V());
    
    // neutral decay to 2gamma, in pi0 rest frame
    if (decay) {
    double gam1_p = CMmomentum (Mass(N0), 0.0, 0.0);
    gamma1.polar(gam1_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 
    gamma1.t(gam1_p);
    gamma2.set(gam1_p,threeVec (0.0, 0.0, 0.0) - gamma1.V());
    gamma2.t(gam1_p);
    }
    /*
     *  compute lorentz factor
     */
    LorentzFactor = (FlatMass || UseBrem) ? 1 : (isobar1_p * cplus_p);
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
	  /* generate vertices */
	  threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;

	// boost gammas to neutral rest frame
	if (decay) {
	tmp.set(n0.t(),zeroVec - n0.V());
	Boost.set(tmp);
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	}
        // boost from c+ c- rest frame to c+ c- n0(X) rest frame

        tmp.set(isobar1.t(),zeroVec - isobar1.V());
        Boost.set (tmp);
        cplus = Boost * cplus;
        cminus = Boost * cminus;

        // boost from 3 particle meson  rest frame to CM
        tmp.set(resonance.t(),zeroVec - resonance.V());
	Boost.set(tmp);

        isobar1 = Boost * isobar1;
        cplus = Boost * cplus;
	cminus = Boost * cminus;
        n0 = Boost * n0;
	if (decay) {
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	}
 

        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        isobar1 = Boost * isobar1;
        cplus = Boost * cplus;
        cminus = Boost * cminus;
        n0 = Boost * n0;
	if (decay) {
	gamma1 = Boost * gamma1;
	gamma2 = Boost * gamma2;
	}



	
	

	if (Print) {

	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "recoil: \n  ";
	  recoil.print();
	  std::cerr << "isobar1\n";
	  std::cerr << "  isobar1 mass: " << isobar1_mass << "\n";
	  std::cerr << "  isobar1 momentum: " << isobar1_p << "\n";
	  std::cerr << "  isobar1:\n    ";
	  isobar1.print();
	  std::cerr << "c+ :\n  ";
	  cplus.print();
	  std::cerr << "n0:\n  ";
	  n0.print();
	  std::cerr << "c-:\n  ";
	  cminus.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << resonance_mass << " " << isobar1_mass*isobar1_mass  << " ";
	  std::cout << ~(n0 + cminus + cplus) << " ";
	  std::cout << pow (~(n0 + cplus), 2.0) << " "
	       << pow (~(n0 + cminus), 2.0) << " "
	       << pow (~(cplus + cminus), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << (n0 + cminus + cplus + recoil).lenSq() << " ";
	  std::cout << (beam - cminus - cplus - gamma1 - gamma2).lenSq() << " ";
	  std::cout << std::endl;
	} 
	else if (debug) {
	  std::cout << "Z " << production.x() << " " << production.y() << " " << production.z() << std::endl;
	}
	else if (txt2part_style){
	  if (decay && printGamma) {
	    std::cout << "5" << std::endl;
	  } 
	  else {
	    std::cout << "4" << std::endl;
	  }
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Baryon, production, recoil);
	  pParticle_txt2part(Cminus, production, cminus);
	  pParticle_txt2part(Cplus, production, cplus);
	  if (decay && printGamma){
	    pParticle_txt2part(Gamma, production, gamma1);
	    pParticle_txt2part(Gamma, production, gamma2);
	  }
	  else{
	    pParticle_txt2part(N0,production,n0);
	  }
	}
	else if (gamp) {

	event evt;
	  
	particle gBeam(PDGtable.get(pid2name(Beam)),Q(Beam));
	particle gTarget(PDGtable.get(pid2name(Target)),Q(Target));
	particle gCminus(PDGtable.get(pid2name(Cminus)),Q(Cminus));
	particle gCplus(PDGtable.get(pid2name(Cplus)),Q(Cplus));
	particle gGamma1(PDGtable.get(pid2name(Gamma)),Q(Gamma));
	particle gGamma2(PDGtable.get(pid2name(Gamma)),Q(Gamma));
	particle gN0(PDGtable.get(pid2name(N0)),Q(N0));

	gBeam.set4P(beam);
	gTarget.set4P(target);
	gCminus.set4P(cminus);
	gCplus.set4P(cplus);
	if (decay) {
	gGamma1.set4P(gamma1);
	gGamma2.set4P(gamma2);
	}
	gN0.set4P(n0);

	evt.beam(gBeam);
	evt.target(gTarget);
	evt.addfinal(gCminus);
	evt.addfinal(gCplus);
	if (decay && printGamma) {
	  evt.addfinal(gGamma1);
	  evt.addfinal(gGamma2);
	}
	else if (printAll) {
	  if (decay) {
	  evt.addfinal(gGamma1);
	  evt.addfinal(gGamma2);
	  }
	  evt.addfinal(gN0);
	}
	else {
	  evt.addfinal(gN0);
	}
	  
	std::cout << evt;
	}
      else {	  	  
	/*
	 *  write event
	 */
	

	if (printBeam)
	  pParticle(Beam,production,beam);


	if (printBaryon)
	  pParticle(Proton,production,recoil);
	pParticle(Cplus,production,cminus);
	pParticle(Cminus,production,cplus);
	if (decay && printGamma) {
	  pParticle(Gamma,production,gamma1);
	  pParticle(Gamma,production,gamma2);
	}
	else {
	  pParticle(N0,production,n0);
	}
      }
      nevents++;
      if (verbose) {
	if (!(nevents % 100)) 
	  std::cerr << nevents << "\r" << flush;
      }
      maxevents--;
      }

    
    }
  }

  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

void bcpcm(int argc, char *argv[],Particle_t Beam,Particle_t Part1,Particle_t Part2,Particle_t Part3)
{
  int icount = 0;
  int Print = 0,
    maxevents = 9999999,
    nevents = 0,
    lfevents = 5000;
  double 
    masslow = Mass(Part2) + Mass(Part3),
    masshigh = 0.,
    t_max,
    expt_max,
    LorentzFactor = 0,
    lfmax = 0,
    resonance_mass,
    isobar1_mass;
  fourVec 
    beam,
    target,
    resonance,
    leading,
    part2,
    part3;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  vector3_t vbeam,pbeam;
  float beamMass;
  int printBaryon = 0;
  int debug = 0;  

  Particle_t Target = Proton;

  /* generate vertices */
  threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
      case 'm':
	ptr++;
	maxevents = atoi (ptr);
	std::cerr << "maxevents: " << maxevents << "\n";
	break;
      case 'l':
	ptr++;
	lfevents = atoi (ptr);
	std::cerr << "lfevents: " << lfevents << "\n";
	break;
      case 'L':
	ptr++;
	masslow = atof(ptr);
	break;
      case 'U':
	ptr++;
	masshigh = atof(ptr);
	break;
      case 'p':
	Print = 1;
	break;
      case 'M':
	break;
      case 'b':
	beamMass = atof(++ptr);
	break;
      case 'B':
	printBaryon = 1;
	break;
      case 'D':
	debug = 1;
	break;
      case 'h':
	Usagecpcm(argv[0]);
	return;
      default:
	std::cerr << "unrecognized argument: " << *ptr << std::endl;
	break;
      }
    }
  }

  double slope = (Slope > 0.0) ? Slope : 4.0;
  int qtot = Q(Beam) + Q(Proton) - Q(Part1) - Q(Part2) - Q(Part3);
  double Threshold =  Mass(Part2) + Mass(Part3);

   beamMass = Mass(Beam);
 
  if(masslow < Threshold)
    masslow = Threshold;

 

  while (maxevents) {
    /*
     *-- use real beam distribution
     */
    generateBeamMCinsideTarget(&vbeam, &pbeam);

    /*
     *-- beam and target in lab frame
     */


    beam = fourVec(sqrt(pow((double) pbeam.x,2.0) + pow((double) pbeam.y,2.0) + pow((double) pbeam.z,2.0) + pow((double) beamMass, 2.0)),threeVec(pbeam.x, pbeam.y, pbeam.z));
    target = fourVec (TARGET_MASS,threeVec(0.0,0.0,0.0));

    /*
     *-- put them into the center of mass frame
     */
    Boost.set (beam + target);
    fourVec CMbeam = Boost * beam;
    fourVec CMtarget = Boost * target;
    double CMenergy = (CMbeam + CMtarget).t();
    double cmphi;

    /*
     *-- generate the resonance and isobar
     */

    /* in case it was not set */
    if(masshigh == 0.)
      masshigh = CMenergy - Mass(Part1);
    if (masshigh < Threshold) {
      std::cerr << "High mass below 2 particle threshold" << std::endl;
      exit(1);
    }

    do {
      resonance_mass = randm(masslow, masshigh);
    } while (resonance_mass < Threshold);

    double beam_p      = CMmomentum (CMenergy, Mass(Beam), PROTON_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, Mass(Part1));
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double costhetax,tx;
    double t;
    double mTot = beamMass + TARGET_MASS + Mass(Part1) + Mass(Part2) - Mass(Part3);
    double s =  pow(sqrt(beam_p * beam_p + beamMass * beamMass)  + sqrt(beam_p * beam_p + TARGET_MASS *TARGET_MASS),2);

    tMin = tmin(pbeam.z,Mass(Beam),PROTON_MASS,Mass(Part1),resonance_mass);

    if (Isotropic)
      costheta = randm(-1.0,1.0);


    else {

      t = getT(tMin,slope);

      //  t = tMin;  // debug
      costheta = 1 + (t - tMin)/(2 * resonance_p * beam_p);
      tx = pow(Mass(Beam),2) + resonance_mass * resonance_mass - 2.0 * resonance_E * sqrt(pow(Mass(Beam),2) + beam_p * beam_p) + 2.0 * resonance_p * beam_p;
  
      costhetax = (t - pow(Mass(Beam),2) - pow(resonance_mass,2) + 2.0 * sqrt(pow(Mass(Beam),2) + beam_p * beam_p) * resonance_E)/(2.0 * resonance_p * beam_p);

    }
    
    
    leading.polar(resonance_p,acos(costheta),cmphi = randm(-M_PI,M_PI));
    leading.t(sqrt(resonance_p * resonance_p + Mass(Part1) * Mass(Part1)));



    /*
     *-- recoil particle
     */
    resonance.set(sqrt(leading.V().lenSq() +  pow (resonance_mass, 2.0)),zeroVec - leading.V());


    if (Print) {
      std::cerr  << "Resonance CM ";
      resonance.print();
      std::cerr << "Leading CM ";
      leading.print();
    }



    /*
     *  now do decay in isobar1 rest frame
     */
    // c+ c-
    double cplus_p = CMmomentum (resonance_mass, Mass(Part2),Mass(Part3));

    part2.polar( cplus_p,acos (randm (-0.999999, 0.999999)),randm (-M_PI, M_PI)); 

    part2.t(sqrt (part2.V().lenSq () + pow (Mass(Part2), 2.0)));
    part3.set(sqrt (part2.V().lenSq () + pow (Mass(Part3), 2.0)),zeroVec - part2.V());
    

    /*
     *  compute lorentz factor
     */
    LorentzFactor = FlatMass ? 1 : resonance_p * cplus_p;
    if (lfevents-- > 0)
      lfmax = LorentzFactor > lfmax ? LorentzFactor : lfmax;
    else {
      if (LorentzFactor > randm (0.0, lfmax)) {
	/* transform all 4-vectors back to lab frame */
        fourVec tmp;


        // boost from c+ c- rest frame to CM rest frame

        tmp.set(resonance.t(),zeroVec - resonance.V());
        Boost.set (tmp);
        part2 = Boost * part2;
        part3 = Boost * part3;


        // boost from CM to target rest frame (lab)
        Boost.set (CMtarget);
        resonance = Boost * resonance;
        leading = Boost * leading;
        part2 = Boost * part2;
        part3 = Boost * part3;

        /* generate vertices */
        threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);
	
	

	if (Print) {
	  std::cerr << "\n\n*** New Event\n";
	  std::cerr << "Beam:\n  ";
	  beam.print();
	  std::cerr << "Beam in CM:\n  ";
	  CMbeam.print();
	  std::cerr << "Target in CM:\n  " ;
	  CMtarget.print();
	  std::cerr << "Resonance\n";
	  std::cerr << "  Resonance mass: " << resonance_mass << "\n";
	  std::cerr << "  Resonance CMmomentum: " << resonance_p << "\n";
	  std::cerr << "  t: " << t << "\n";
	  std::cerr << "  Resonance:\n "    ;
	  resonance.print();
	  std::cerr << "leading: \n  ";
	  leading.print();
	  std::cerr << "c+ :\n  ";
	  part2.print();
	  std::cerr << "c-:\n  ";
	  part3.print();
	  //	   std::cerr << "vertices:\n";
	  //	   std::cerr << "  prod: " << production;
	  std::cerr << "Lorentz factor: " << LorentzFactor << "\n";
	  std::cerr << "icount: " << icount << std::endl;
	}
	// calculate masses for dalitz plots
	if (dalitz) {
	  std::cout << "DALITZ ";
	  std::cout << resonance_mass << " ";
	  std::cout << ~(leading + part3 + part2) << " ";
	  std::cout << pow (~(leading + part2), 2.0) << " "
	       << pow (~(leading + part3), 2.0) << " "
	       << pow (~(part2 + part3), 2.0) << " ";
	  std::cout << masslow << " " << masshigh << " ";
	  std::cout << std::endl;
	} 
	else if (debug) {
	  std::cout << "T " << (beam - part2 - part3).lenSq() << " " << (target - leading).lenSq() << " " << t << " " << tMin << " " << costheta << " " <<  costhetax << " " << " " << cmphi << " " << (t - tMin) << " " << tx << " " << std::endl;
	}
	if (txt2part_style){
	  std::cout << "3" << std::endl;
	  std::cout << EbeamZ << " " << beam.t() << " " << -production.z()/SPEED_OF_LIGHT << std::endl;; 
	  pParticle_txt2part(Part1, production, leading);
	  pParticle_txt2part(Part3, production, part3);
	  pParticle_txt2part(Part2, production, part2);
	}
	else if (gamp) {
	  event evt;
	  
	  particle gBeam(PDGtable.get(pid2name(Beam)),Q(Beam));
	  particle gTarget(PDGtable.get(pid2name(Target)),Q(Target));
	  particle gLeading(PDGtable.get(pid2name(Part1)),Q(Part1));
	  particle gPart3(PDGtable.get(pid2name(Part3)),Q(Part3));
	  particle gPart2(PDGtable.get(pid2name(Part2)),Q(Part2));

	  gBeam.set4P(beam);
	  gTarget.set4P(target);
	  gLeading.set4P(leading);
	  gPart3.set4P(part3);
	  gPart2.set4P(part2);

	  evt.beam(gBeam);
	  evt.target(gTarget);
	  evt.addfinal(gLeading);
	  evt.addfinal(gPart3);
	  evt.addfinal(gPart2);
	  
	  std::cout << evt;
	  
	  
	}
	else {
	  /*
	   *  write event
	   */
	  nevents++;
	  if (verbose) {
	    if (!(nevents % 100)) 
	      std::cerr << nevents << "\r" << flush;
	  }

	  if (printBeam)
	    pParticle(Beam,production,beam);


	  pParticle(Part1,production,leading);
	  pParticle(Part2,production,part2);
	  pParticle(Part3,production,part3);

	}

	maxevents--;
      }
    }
  }
  if (verbose)
    std::cerr << nevents << " Events generated" << std::endl;
}

/*----------------End of charged1 charged2 -----------------------------------*/





#include <cpcmEXP.h>

#include <bpn.h>



/*
 * 
 */

