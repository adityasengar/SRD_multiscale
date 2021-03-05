#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define NPART_MAX 1000000    // Maximum number of SRD particles
#define NCOL_MAX  1000       // Maximum number of colloidal particles (not currently used in Code2)
#define LBOX_MAX  100        // Maximum number of collision cells in box size (for any dimension)

// Variables for velocity profile calculations
extern int ybin;             // Number of bins in y-direction for velocity profiles
extern int zbin;             // Number of bins in z-direction for velocity profiles
extern double *vel_sum;     // Array to sum up z-velocities for a particular bin
extern int *counter;         // Integer counter to update by 1 when a particle lies in a bin
extern double dbin;         // Bin size in y-direction
extern double dbin_z;       // Bin size in z-direction

// Variables for correlation functions (MSD, ACF)
extern int numcor;           // Number of correlation steps
extern int num_dtau;         // Time interval between correlation measurements

#endif /* CONSTANTS_H_ */
\n\/\* Fluid density \*\/\n#define RHO 1.0
