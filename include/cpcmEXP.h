/* 
 * File:   cpcmEXP.h
 * Author: dennisweygand
 *
 * Created on March 28, 2012, 8:20 AM
 */

#ifndef CPCMEXP_H
#define	CPCMEXP_H


void cpcmEXP(int argc, char *argv[],Particle_t Beam,Particle_t Cplus,Particle_t Cminus,Particle_t Recoil)
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
  double alpha = 1.0;	

  Particle_t Target = Proton;

  /* generate vertices */
  threeVec production = threeVec (vbeam.x, vbeam.y, vbeam.z);

  float tMin;

  for (int iarg = 1; iarg < argc; ++iarg) {
    char *ptr = argv[iarg];
    if (*ptr == '-') {
      ptr++;
      switch (*ptr) {
	case 'a':
	ptr++;
	alpha = atof(ptr);
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
      case 'D':
	debug = 1;
	break;
      case 'h':
	std::cerr << "\t-a#\t#: slope of mass distribution" << std::endl;
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
      resonance_mass = expDist(alpha);
	if (debug)
	cout << "D " << resonance_mass << endl;
    } while (resonance_mass < masslow ||  resonance_mass > masshigh);

    double beam_p      = CMmomentum (CMenergy, Mass(Beam), PROTON_MASS);
    double resonance_p = CMmomentum (CMenergy, resonance_mass, Mass(Baryon));
    double resonance_E = sqrt (pow (resonance_mass, 2.0) + pow (resonance_p, 2.0));
    double costheta;
    double costhetax,tx;
    double t;
    double mTot = beamMass + TARGET_MASS + Mass(Baryon) + Mass(Cplus) - Mass(Cminus);
    double s =  pow(sqrt(beam_p * beam_p + beamMass * beamMass)  + sqrt(beam_p * beam_p + TARGET_MASS *TARGET_MASS),2);

    tMin = tmin(pbeam.z,Mass(Beam),PROTON_MASS,resonance_mass,Mass(Baryon));

    if (Isotropic)
      costheta = randm(-1.0,1.0);


    else {

      t = getT(tMin,slope);

      //  t = tMin;  // debug
      costheta = 1 + (t - tMin)/(2 * resonance_p * beam_p);
      tx = pow(Mass(Beam),2) + resonance_mass * resonance_mass - 2.0 * resonance_E * sqrt(pow(Mass(Beam),2) + beam_p * beam_p) + 2.0 * resonance_p * beam_p;
  
      costhetax = (t - pow(Mass(Beam),2) - pow(resonance_mass,2) + 2.0 * sqrt(pow(Mass(Beam),2) + beam_p * beam_p) * resonance_E)/(2.0 * resonance_p * beam_p);

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

#endif	/* CPCMEXP_H */

