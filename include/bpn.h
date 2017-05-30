/* 
 * File:   bpn.h
 * Author: dennisweygand
 *
 * Created on March 29, 2012, 7:44 PM
 */

#ifndef BPN_H
#define	BPN_H

void bpn(int argc, char *argv[],Particle_t Beam,Particle_t Part1,Particle_t Part2,Particle_t Part3,Particle_t D1,Particle_t D2)
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
    part3,
    d1,d2;
  lorentzTransform Boost;
  
  threeVec zeroVec = threeVec (0.0, 0.0, 0.0);
  fourVec tmp;
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

    // now decay leading

    double d_p = CMmomentum(Mass(Part1),Mass(D1),Mass(D2));
    double cth = randm(-1,1);
    d1.polar(d_p,acos(cth),randm(-M_PI,M_PI));
    d1.t(sqrt(d_p*d_p +Mass(D1) * Mass(D1)));

    d2.set(sqrt(d_p * d_p + Mass(D2) * Mass(D2)),zeroVec - d1.V());

    // transform
    tmp.set(leading.t(),zeroVec - leading.V());
    Boost.set(tmp);
    d1 *= Boost;
    d2 *= Boost;



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
	d1 *= Boost;
	d2 *= Boost;
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
	  particle gD1(PDGtable.get(pid2name(D1)),Q(D1));
	  particle gD2(PDGtable.get(pid2name(D2)),Q(D2));
          
          cout << pid2name(Part2) << " " << pid2name(Part3) <<  endl;

	  gBeam.set4P(beam);
	  gTarget.set4P(target);
	  gLeading.set4P(leading);
	  gPart3.set4P(part3);
	  gPart2.set4P(part2);
	  gD1.set4P(d1);
	  gD2.set4P(d2);

	  evt.beam(gBeam);
	  evt.target(gTarget);
//	  evt.addfinal(gLeading);
	  evt.addfinal(gPart3);
	  evt.addfinal(gPart2);
	  evt.addfinal(gD1);
	  evt.addfinal(gD2);
	  
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


#endif	/* BPN_H */

