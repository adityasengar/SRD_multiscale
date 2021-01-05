#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// Include custom header files
#include "constants.h"
#include "variables.h"
#include "initialise.h"
#include "functions.h"
#include "memory.h" // Include memory.h for memory management functions

// Global variables (declared in variables.h, defined here or in other .c files)
// Simulation parameters
int Npart, Nsrd, NA;
double *Rx, *Ry, *Rz;
double *Vx, *Vy, *Vz;
double *Rz_artificial;

double *Rx_adsA, *Ry_adsA, *Rz_adsA;
double *Rx_adsB, *Ry_adsB, *Rz_adsB;

int *list, *list_mass;
int L_x, L_y, L_z;
int ***head, ***head_mass;
int ***nearwalls;
int doshift;

FILE *fp_traj;
double box_x, box_y, box_z;
double t, dt, tmin, tmax, tstart, tau;
int dt_flow;
double shiftx, shifty, shiftz;
double kTmeas, kTo;
double conc;
int walls;
int therm;
int CONV_DIFF;
int CONV_DIFF_mass;
int x1, x2;
double alpha, beta;
int gamma1;
double t_step, t_interval;
double g;
int open; // 1: open system in z-direction, 0: periodic boundaries in z-direction

// MSD and ACF variables
double *cor;
double *cor1;
double *cor2;
int *ACFcount;
int *ACFcount1;
int *ACFcount2;
double **storex;
double **storey;
double **storez;
double **storex1;
double **storey1;
double **storez1;
double **storex2;
double **storey2;
double **storez2;

// Slit tracking variables
int track;
int *tracker;
int num_track;

// Velocity profile variables
double *velz_meas;
double *velz_flow;
double *N_average;
int velz_bin;
int ybin, zbin;
double *vel_sum;
int *counter;
double dbin, dbin_z;
int numcor, num_dtau;
double *y_sum;
double *z_sum;
double *z_sum1;
int N_recur;

// Adsorption-desorption reaction model variables
int adsA, adsB;
int catalytic_sites;
double p_adsA, p_adsB;
double p_desA, p_desB;
double p_react;
double thetaA, thetaB, theta;
int **catalyst;
double timey;

// For homogenous function
int *arrays;


void initial_conditions() {
    time_t t_seed;
    srand((unsigned) time(&t_seed)); // Seed random number generator
    int i;
    double sum_vx = 0.0, sum_vy = 0.0, sum_vz = 0.0;
    double total_mass = 0.0;

    // Initialize SRD particles
    if (CONV_DIFF == 1) { // Initial rectangle concentration method
        for (i = 0; i < Nsrd; i++) {
            double r3_rand = box_z * (double)rand() / (double)RAND_MAX;
            double r1_rand = box_x * (double)rand() / (double)RAND_MAX;
            double r2_rand;

            if (r3_rand < alpha * box_z) {
                double rand_i = (double)rand() / (double)RAND_MAX;
                if (rand_i < (double)x1 * beta / (double)gamma1) {
                    r2_rand = beta * box_y * (double)rand() / (double)RAND_MAX;
                } else if (rand_i < ((double)x1 * beta / (double)gamma1 + (1.0 - 2.0 * beta) * (double)x2 / (double)gamma1)) {
                    r2_rand = beta * box_y + (1.0 - 2.0 * beta) * box_y * (double)rand() / (double)RAND_MAX;
                }
                else {
                    r2_rand = (1.0 - beta) * box_y + beta * box_y * (double)rand() / (double)RAND_MAX;
                }
            } else {
                r2_rand = box_y * (double)rand() / (double)RAND_MAX;
            }

            Rx[i] = r1_rand;
            Ry[i] = r2_rand;
            Rz[i] = r3_rand;
            Rz_artificial[i] = r3_rand;

            Vx[i] = sqrt(kTo) * neargauss();
            Vy[i] = sqrt(kTo) * neargauss();
            Vz[i] = sqrt(kTo) * neargauss();

            sum_vx += Vx[i];
            sum_vy += Vy[i];
            sum_vz += Vz[i];
            total_mass += 1.0;
        }
        // Initialize foreign particles (A)
        for (i = Nsrd; i < Npart; i++) {
            double r1_rand = box_x * (double)rand() / (double)RAND_MAX;
            double r2_rand = box_y * (double)rand() / (double)RAND_MAX;
            double r3_rand = box_z * (double)rand() / (double)RAND_MAX;

            if (r3_rand < alpha * box_z) {
                double rand1 = (double)rand() / (double)RAND_MAX;
                if (rand1 >= 0.5) {
                    r2_rand = beta * box_y * (double)rand() / (double)RAND_MAX;
                } else {
                    r2_rand = (1.0 - beta) * box_y + beta * box_y * (double)rand() / (double)RAND_MAX;
                }
            } else {
                r2_rand = box_y * (double)rand() / (double)RAND_MAX;
            }

            Rx[i] = r1_rand;
            Ry[i] = r2_rand;
            Rz[i] = r3_rand;

            Vx[i] = sqrt(kTo / massA) * neargauss();
            Vy[i] = sqrt(kTo / massA) * neargauss();
            Vz[i] = sqrt(kTo / massA) * neargauss();

            sum_vx += massA * Vx[i];
            sum_vy += massA * Vy[i];
            sum_vz += massA * Vz[i];
            total_mass += massA;
        }
    } else { // Default initialization (homogeneous distribution)
        for (i = 0; i < Nsrd; i++) {
            Rx[i] = box_x * (double)rand() / (double)RAND_MAX;
            Ry[i] = box_y * (double)rand() / (double)RAND_MAX;
            Rz[i] = box_z * (double)rand() / (double)RAND_MAX;
            Rz_artificial[i] = Rz[i]; // Initialize artificial z-position

            Vx[i] = sqrt(kTo) * neargauss();
            Vy[i] = sqrt(kTo) * neargauss();
            Vz[i] = sqrt(kTo) * neargauss();

            sum_vx += Vx[i];
            sum_vy += Vy[i];
            sum_vz += Vz[i];
            total_mass += 1.0;
        }
        // Initialize foreign particles (A)
        for (i = Nsrd; i < Npart; i++) {
            Rx[i] = box_x * (double)rand() / (double)RAND_MAX;
            Ry[i] = box_y * (double)rand() / (double)RAND_MAX;
            Rz[i] = box_z * (double)rand() / (double)RAND_MAX;

            Vx[i] = sqrt(kTo / massA) * neargauss();
            Vy[i] = sqrt(kTo / massA) * neargauss();
            Vz[i] = sqrt(kTo / massA) * neargauss();

            sum_vx += massA * Vx[i];
            sum_vy += massA * Vy[i];
            sum_vz += massA * Vz[i];
            total_mass += massA;
        }
    }

    // Remove center of mass motion
    sum_vx /= total_mass;
    sum_vy /= total_mass;
    sum_vz /= total_mass;

    for (i = 0; i < Npart; i++) {
        if (i < Nsrd) { // SRD particles have mass 1.0
            Vx[i] -= sum_vx;
            Vy[i] -= sum_vy;
            Vz[i] -= sum_vz;
        } else { // Foreign particles have mass massA
            Vx[i] -= sum_vx / massA;
            Vy[i] -= sum_vy / massA;
            Vz[i] -= sum_vz / massA;
        }
    }
    COM_vel_zero(); // Further adjust COM velocity if needed (implementation in calc.c)
}


int main( int argn, char * argv[] )
{
    // Initialize adsorption counters and probabilities
    adsA = 0;
    adsB = 0;
    p_adsA = 0.5;
    p_adsB = 0.0;
    p_react = 0.0;
    p_desA = 0.5;
    p_desB = 0.0;
    open = 1; // 1: open system in z-direction, 0: periodic boundaries in z-direction

    // Open trajectory file
    fp_traj = fopen("trajectories.pdb", "w");
    if (fp_traj == NULL) {
        fprintf(stderr, "Error: Could not open trajectories.pdb for writing.\n");
        return EXIT_FAILURE;
    }

    // Simulation setup parameters
    gamma1 = 5;  // Mass density per unit cell (or similar parameter)
    int which_particle_count_method = 0; // 0: mass density is constant per unit cell, 1: particles constant
    track = 0;  // 1: track slit behavior in 3D non-boundary medium, 0: off
    walls = 1; // 1: solid walls present, 0: periodic boundaries
    therm = 1; // 1: turn on Galilean thermostat, 0: off
    g = 0.000; // External force (gravity/acceleration)

    // Box dimensions (number of collision cells)
    L_x = 5;
    L_y = 5;
    L_z = 10;

    // Parameters for specific simulation methods
    alpha = 0.02; // Fraction of L_z used to insert particles for axial dispersion method
    beta = 0.2;   // Parameter for slit geometry (e.g., width of slit)
    x1 = 22;      // Parameters for different methods (not used in current axial dispersion)
    x2 = 2;       // Parameters for different methods (not used in current axial dispersion)

    // Simulation modes for convection-diffusion
    CONV_DIFF = 0; // 0: MSD/ACF, 1: Initial rectangle concentration, 2: Axial dispersion method
    CONV_DIFF_mass = 0; // 0: Same mass particles (axial dispersion), 1: Different mass particles

    int Neff = L_x * L_y * L_z * gamma1; // Effective number of particles

    // Diffusion calculation method
    // 0: MSD technique, 1: Velocity Autocorrelation (dual), 2: Stationary velocity method,
    // 3: Unary velocity autocorrelation, 4: Green-Kubo (something function)
    int diffusion_method = 0;

    // Loop for multiple simulation runs (if needed, currently set to run once)
    for (int j = 0; j < 1; j++) // Original loop was j<5, j=j+20, which runs once.
    {
        massA = 1.0; // Mass of foreign particles (A)
        conc = 0.0;  // Concentration of foreign particles (A)

        for (conc = 0.0; conc <= 0.1; conc += 1.1) // Original loop was conc<=1.0, conc+=1.1, which runs once.
        {
            dt = 0.1;   // Time step for simulation
            thetaA = 0.0; // Initial surface coverage of A
            thetaB = 0.0; // Initial surface coverage of B
            theta = 1.0 - thetaA - thetaB; // Initial vacant sites

            // Calculate number of catalytic sites (example calculation, adjust as needed)
            catalytic_sites = Npart * 0.2; // Different calculation than Code1
            printf("Number of catalytic sites: %d\n", catalytic_sites);

            dt_flow = (int)(10 * dt); // Time interval for flow updates

            // Determine total number of particles based on chosen method
            if (which_particle_count_method == 0) { // Mass density is constant per unit cell
                Npart = (int)((double)Neff / (conc * massA + 1 - conc));
            } else {
                Npart = Neff;
            }
            Nsrd = (int)((1 - conc) * Npart); // Number of SRD particles
            NA = Npart - Nsrd; // Number of foreign particles (A)

            tmin = 1000.0 * dt;  // Start of evaluation for various parameters
            tmax = 100.0 * dt * 1000.0; // End of simulation time
            t_step = 1.0 * dt;   // Time pulse width for axial dispersion
            t_interval = 1000.0 * dt; // Time between which next pulse is sent
            tstart = 0.0;        // Simulation start time

            int trajoutput_interval = 50; // Interval for recording trajectories
            int print_energy_interval = 1000; // Interval for printing energy/momentum

            numcor = 201; // Number of correlation steps for ACF/MSD
            num_dtau = 1; // Time interval between correlation measurements

            allocateMemory(); // Allocate memory for all arrays
            if (CONV_DIFF_mass == 1 && walls == 1) {
                alloc_memory_CONV_DIFF_mass(Npart);
            }

            kTo = 1.0; // Target thermal energy
            tau = 1.0; // Scaling factor (from calc.c line 339,340 comment)

            int slit_trajout_interval = 10; // Interval for slit tracking trajectory output

            box_x = (double)L_x;
            box_y = (double)L_y;
            box_z = (double)L_z;

            printf("\nNpart: %d, tmax: %f\n", Npart, tmax);
            t = tstart;
            int tstep = 0;
            kTmeas = 1.0; // Initial measured thermal energy

            initial_conditions(); // Set up initial particle positions and velocities

            alloc_memory_ACF(Npart, numcor); // Allocate memory for ACF/MSD arrays

            // Initial calls for data collection functions based on simulation mode
            if (walls == 0 && track == 0) {
                if (diffusion_method == 0) updatecorrelator_MSD(1, massA);
                if (diffusion_method == 1) dual_velocity_correlator(1, massA);
                if (diffusion_method == 2) new_diff(1, t);
                if (diffusion_method == 3) velocity_correlator(1, dt);
                if (diffusion_method == 4) something(1, massA);
            } else if (walls == 1 && track == 0 && CONV_DIFF == 0) {
                velocityprofile(0, t);
            } else if (walls == 1 && CONV_DIFF == 1) {
                profile_maker(0, t);
                update_velz_bin(0, t);
                velocityprofile(0, t);
            }

            if (walls == 1 && CONV_DIFF == 2) {
                homogenous(0, (double)((int)t / dt) % (int)(t_interval / dt));
                velocityprofile(0, t);
            }

            if (track == 1) slit_tracker(0, t); // Pass t to slit_tracker

            // Determine which cells overlap with walls
            if (walls == 1) determinenearwalls();

            // Record initial positions if interval is set
            if (trajoutput_interval > 0) {
                // record_trajectories(); // This was commented out in original, keeping it commented.
            }

            int total_timesteps = (int)(tmax / dt);

            // MAIN SIMULATION LOOP OVER TIME
            while (t < tmax) {
                // Progress indicators
                if (tstep == total_timesteps / 10) printf("10%....\n");
                if (tstep == 2 * total_timesteps / 10) printf("20%....\n");
                if (tstep == 3 * total_timesteps / 10) printf("30%....\n");
                if (tstep == 4 * total_timesteps / 10) printf("40%....\n");
                if (tstep == 5 * total_timesteps / 10) printf("50%....\n");
                if (tstep == 6 * total_timesteps / 10) printf("60%....\n");
                if (tstep == 7 * total_timesteps / 10) printf("70%....\n");
                if (tstep == 8 * total_timesteps / 10) printf("80%....\n");
                if (tstep == 9 * total_timesteps / 10) printf("90%....\n");

                // Record trajectories
                if (walls == 1 && tstep % trajoutput_interval == 0 && t > 0 && track == 0 && CONV_DIFF == 0) {
                    record_trajectories();
                }

                // Update flow profiles for CONV_DIFF == 1
                if (walls == 1 && CONV_DIFF == 1 && tstep % dt_flow == 0 && t > tmin) {
                    profile_maker(1, t);
                    update_velz_bin(1, t);
                    velocityprofile(1, t);
                }

                doshift = 1; // Enable random grid shift

                updatepositions(); // Update particle positions
                adsorption();      // Handle adsorption/desorption/reaction on catalytic surface
                createlinkedlist(); // Create list of particles in (possibly shifted) cells
                collide_particles(); // Perform particle collisions
                // COM_vel_zero(); // Commented out in original, keeping it commented.
                updatevelocities(); // Update velocities due to external fields

                // Slit tracking
                if (walls == 0 && track == 1) {
                    if (tstep % slit_trajout_interval == 0) {
                        slit_tracker(1, t);
                    }
                }

                // Optionally output energy and momentum
                // if (print_energy_interval > 0 && tstep % print_energy_interval == 0) {
                //     energy_momentum(0, t);
                // }

                // Update correlators for MSD/ACF
                if (walls == 0 && track == 0) {
                    if (tstep % num_dtau == 0 && t > tmin) {
                        if (diffusion_method == 0) updatecorrelator_MSD(0, massA);
                        if (diffusion_method == 1) dual_velocity_correlator(0, massA);
                        if (diffusion_method == 2) new_diff(0, t);
                        if (diffusion_method == 3) velocity_correlator(0, dt);
                        if (diffusion_method == 4) something(0, massA);
                    }
                }

                // Homogeneous updates for CONV_DIFF == 2
                if (walls == 1 && CONV_DIFF == 2 && t >= tmin) {
                    homogenous(1, (double)((int)t / dt) % (int)(t_interval / dt));
                }

                tstep++;
                t = t + dt;
            }

            // Final calls for data collection functions after simulation loop
            velocityprofile(2, t);
            homogenous(2, (double)((int)t / dt) % (int)(t_interval / dt));
            printf("Finalizing data output...\n");

            if (fp_traj != NULL) {
                fclose(fp_traj);
            }

            if (walls == 1 && CONV_DIFF == 0) {
                velocityprofile(2, t);
            } else if (CONV_DIFF == 1) {
                profile_maker(2, t);
            } else {
                if (diffusion_method == 0) updatecorrelator_MSD(2, massA);
                if (diffusion_method == 1) dual_velocity_correlator(2, massA);
                if (diffusion_method == 2) new_diff(2, t);
                if (diffusion_method == 3) velocity_correlator(2, dt);
                if (diffusion_method == 4) something(2, massA);
            }

            // Free memory
            deallocateMemory();
            if (walls == 0) { // Only deallocate ACF memory if it was allocated (walls==0 implies MSD/ACF)
                deallocateMemory_ACF();
            }
        }
    }
    printf("Simulation finished!\n");
    return 0;
}
