#ifndef VARIABLES_H_
#define VARIABLES_H_

#include <stdio.h>
#include "constants.h" // For NPART_MAX, LBOX_MAX

// Simulation parameters
extern int Npart, Nsrd, NA;            // Number of total particles, SRD particles, and foreign (A) particles
extern double *Rx, *Ry, *Rz;           // SRD particle positions
extern double *Vx, *Vy, *Vz;           // SRD particle velocities
extern double *Rz_artificial;          // Artificial z-position for axial dispersion method

// Adsorbed particle positions (for catalytic reactions)
extern double *Rx_adsA, *Ry_adsA, *Rz_adsA;
extern double *Rx_adsB, *Ry_adsB, *Rz_adsB;

// Linked-list for collision cells
extern int *list;                      // Linked-list array for SRD particles
extern int *list_mass;                 // Linked-list array for particles with different masses (Axial Dispersion)
extern int ***head;                    // Head array of linked-list for SRD particles
extern int ***head_mass;               // Head array of linked-list for particles with different masses

// Simulation box dimensions and grid
extern int L_x, L_y, L_z;              // Number of collision cells in x, y, z dimensions
extern double box_x, box_y, box_z;     // Real box size in x, y, z dimensions

// Time parameters
extern double t, dt, tmin, tmax, tstart, tau; // Current time, time step, min/max simulation time, start time, scaling factor
extern int dt_flow;                    // Time grid for flow-related updates
extern double t_step;                 // Time pulse width for axial dispersion
extern double t_interval;             // Time between which next pulse is sent

// Random grid shift
extern int doshift;                    // Flag for random grid shift (1: apply, 0: no shift)
extern double shiftx, shifty, shiftz;  // Random grid shift values [0,1]

// Thermostat and external forces
extern double kTmeas, kTo;             // Measured and target thermal energy
extern int therm;                      // Flag to turn on Galilean thermostat (1: on, 0: off)
extern double g;                      // External force (gravity/acceleration)

// Boundary conditions
extern int walls;                      // Flag for boundary type (1: solid walls, 0: periodic boundaries)
extern int ***nearwalls;               // 1 if cell is near a wall
extern int open;                       // 1: open system in z-direction, 0: periodic boundaries in z-direction

// Concentration and mass parameters
extern double conc;                   // Concentration of foreign particles (A)
extern double massA;                  // Mass of foreign particles (A)

// Simulation modes/methods
extern int CONV_DIFF;                  // 0: MSD/ACF, 1: Initial rectangle concentration, 2: Axial dispersion method
extern int CONV_DIFF_mass;             // 0: Same mass particles (axial dispersion), 1: Different mass particles

// Parameters for specific methods (x1, x2, alpha, beta, gamma1)
extern int x1, x2;                     // Used for different methods (not used in current axial dispersion)
extern double alpha;                  // Fraction of L_z used to insert particles for axial dispersion
extern double beta;                   // Parameter for slit geometry
extern int gamma1;                     // Mass density per unit cell (or similar parameter)

// File pointers
extern FILE *fp_traj;                  // File pointer for trajectory output

// MSD and ACF variables
extern double *cor;
extern double *cor1;
extern double *cor2;
extern int *ACFcount;
extern int *ACFcount1;
extern int *ACFcount2;
extern double **storex;
extern double **storey;
extern double **storez;
extern double **storex1;
extern double **storey1;
extern double **storez1;
extern double **storex2;
extern double **storey2;
extern double **storez2;

// Slit tracking variables
extern int track;                      // Flag to track behavior of a slit (1: on, 0: off)
extern int *tracker;                   // Array to store tracked particle indices
extern int num_track;                  // Number of tracked particles

// Velocity profile variables (for update_velz_bin and velocity_changer)
extern double *velz_meas;              // Measured z-velocities in bins
extern double *velz_flow;              // Flow z-velocities in bins
extern double *N_average;              // Number of particles averaged in each bin
extern int velz_bin;                  // Number of bins for z-velocity profile
extern double N_recur;                 // Variable for velocity updation when entering the system again

// Adsorption-desorption reaction model variables
extern int adsA, adsB;                 // Number of adsorbed A and B particles
extern int catalytic_sites;            // Total number of catalytic sites
extern double p_adsA, p_adsB;         // Probability of adsorption for A and B
extern double p_desA, p_desB;         // Probability of desorption for A and B
extern double p_react;                // Probability of reaction (A -> B)
extern double thetaA, thetaB, theta;  // Surface coverage of A, B, and vacant sites
extern int **catalyst;                 // Array to mark catalyst regions and particle history
extern double timey;                  // Temporary variable for time calculation in updatepositions

// For homogenous function
extern int *arrays;

#endif /* VARIABLES_H_ */
