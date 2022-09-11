/*
 * particleType.h
*/
#pragma once

#define ETA_MASS          0.54745
#define PROTON_MASS       0.93827
#define NEUTRON_MASS      0.93957
#define OMEGA_MASS        0.78194
#define PHI_MASS          1.019413

typedef enum {

    /*
     * These constants are defined to be
     * same as GEANT. see http://wwwcn.cern.ch/asdoc/geant/H2GEANTCONS300.html
     * for more details.
    */

    Gamma = 1,
    Positron = 2,
    Electron = 3,
    Pi0 = 7,
    PiPlus = 8,
    PiMinus = 9,
    KLong = 10,
    KPlus = 11,
    KMinus = 12,
    Neutron = 13,
    Proton = 14,
    AntiProton = 15,
    KShort = 16,
    Eta = 17,
    AntiNeutron = 25,
    Alpha = 47,

    Rho0 = 57,
    omega = 60,
    EtaPrime = 61,
    phiMeson = 62

} Particle_t;


#ifndef BIT
#define BIT(n) (1<<(n))
#endif
