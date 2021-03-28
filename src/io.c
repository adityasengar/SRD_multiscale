#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For strlen, snprintf
#include <math.h>   // For floor, round
#include "initialise.h"
#include "variables.h"
#include "constants.h" // For NUM_PARTICLES_MAX, LBOX_MAX
#include "memory.h" // For memory allocation functions

// Global variables for trajectory recording and correlation functions
static int frame = 0;
static double average = 0.0;
static double average1 = 0.0;
static int cc = 0;
static int cc1 = 0;

/**
 * @brief Records particle trajectories in PDB format to a file.
 * This function is called periodically to save the state of the simulation
 * for visualization in tools like VMD.
 */
void record_trajectories()
{
  int i;

  fprintf(fp_traj, "REMARK TIME = %f\n", t);
  fprintf(fp_traj, "CRYST1 %6.2f %6.2f %6.2f    5.0    10.0    10.0 P 1          1\n", box_x, box_y, box_z);

  // Output foreign particles (type S)
  for (i = Nsrd; i < Npart; i++) {
    fprintf(fp_traj, "HETATM %5i  S 1 UNK     1    %7.4f %7.4f %7.4f   1.0   1.0\n", i, Rx[i], Ry[i], Rz[i]);
  }
  // Output adsorbed A particles (type N)
  for (i = 0; i < adsA; i++) {
    fprintf(fp_traj, "HETATM %5i  N 1 UNK     1    %7.4f %7.4f %7.4f   1.0   1.0\n", i, Rx_adsA[i], Ry_adsA[i], Rz_adsA[i]);
  }
  // Output adsorbed B particles (type O)
  for (i = 0; i < adsB; i++) {
    fprintf(fp_traj, "HETATM %5i  O 1 UNK     1    %7.4f %7.4f %7.4f   1.0   1.0\n", i, Rx_adsB[i], Ry_adsB[i], Rz_adsB[i]);
  }
  // Output placeholder particles (type C) - original code had this, purpose unclear without more context
  for (i = 0; i < 200 - adsA - adsB - NA; i++) { // 'NA' here seems to be a particle count, not the global NA
    fprintf(fp_traj, "HETATM %5i  C 1 UNK     1    %7.4f %7.4f %7.4f   1.0   1.0\n", i, 1.0, 1.0, 1.0);
  }
  fprintf(fp_traj, "END\n");
}

/**
 * @brief Periodically checks and outputs total energy and center of mass momentum.
 * This function is used to verify conservation laws in the simulation.
 * @param info Flag to control output (0: print, 1: return).
 * @param time Current simulation time.
 */
void energy_momentum(int info, double time)
{
  double kinetic = 0.0;
  double sum_vx = 0.0, sum_vy = 0.0, sum_vz = 0.0;
  double total_mass = 0.0;

  for (int i = 0; i < Npart; i++) {
    if (i < Nsrd) { // SRD particles (mass 1.0)
        sum_vx += Vx[i];
        sum_vy += Vy[i];
        sum_vz += Vz[i];
        kinetic += (Vx[i] * Vx[i] + Vy[i] * Vy[i] + Vz[i] * Vz[i]);
        total_mass += 1.0;
    } else { // Foreign particles (mass massA)
        sum_vx += massA * Vx[i];
        sum_vy += massA * Vy[i];
        sum_vz += massA * Vz[i];
        kinetic += massA * (Vx[i] * Vx[i] + Vy[i] * Vy[i] + Vz[i] * Vz[i]);
        total_mass += massA;
    }
  }
  kinetic *= 0.5; // Total kinetic energy

  if (total_mass > 0) {
    sum_vx /= total_mass;
    sum_vy /= total_mass;
    sum_vz /= total_mass;
  }

  if (info == 0) {
    // printf( "%f %f %f %f %f %f %f %f %f\n", time, kinetic+kineticcol+Ucs+Ucc, kinetic, kineticcol, Ucs, Ucc, sumvx, sumvy, sumvz );
    // Original code had commented out output, keeping it commented.
    return;
  }
}

/**
 * @brief Updates and calculates Mean Squared Displacement (MSD) correlator.
 * This function is used to calculate the diffusion coefficient.
 * @param reset 1: initialize, 0: update, 2: finalize and write output.
 * @param massA Mass of foreign particles.
 */
void updatecorrelator_MSD(int reset, double massA)
{
    FILE* output_file;
    char filename[256]; // Buffer for filename

    // Construct filename based on concentration and massA
    snprintf(filename, sizeof(filename), "data/Diffusion_concA=%.3f_m=%.1f.dat", conc, massA);

    if (reset == 1) { // Initialize
        frame = 0;
        output_file = fopen(filename, "w");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        fprintf(output_file, "mass of particle A = %f\n", massA);
        fclose(output_file);

        for (int i = 0; i < Npart; i++) {
            storex[i][frame] = Rx[i];
            storey[i][frame] = Ry[i];
            storez[i][frame] = Rz[i];
        }
        frame = 1;
    } else if (reset == 0) { // Update
        int maxcor_val = numcor - 1;
        if (frame < numcor) maxcor_val = frame;
        int curframe = frame % numcor;
        int prevframe = (curframe - 1 + numcor) % numcor;

        for (int i = 0; i < Npart; i++) {
            double deltax = Rx[i] - storex[i][prevframe];
            double deltay = Ry[i] - storey[i][prevframe];
            double deltaz = Rz[i] - storez[i][prevframe];

            // Unfold coordinates for periodic boundaries
            double xunfold = storex[i][prevframe] + deltax - box_x * round(deltax / box_x);
            double yunfold = storey[i][prevframe] + deltay - box_y * round(deltay / box_y);
            double zunfold = storez[i][prevframe] + deltaz - box_z * round(deltaz / box_z);

            storex[i][curframe] = xunfold;
            storey[i][curframe] = yunfold;
            storez[i][curframe] = zunfold;

            for (int icor = 0; icor <= maxcor_val; icor++) {
                int icorframe = (curframe - icor + numcor) % numcor;
                cor[icor] += (xunfold - storex[i][icorframe]) * (xunfold - storex[i][icorframe]) + 
                             (yunfold - storey[i][icorframe]) * (yunfold - storey[i][icorframe]) + 
                             (zunfold - storez[i][icorframe]) * (zunfold - storez[i][icorframe]);
                ACFcount[icor] += 1;
            }
        }
        frame++;
    } else { // Finalize and write output
        output_file = fopen(filename, "a");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        for (int icor = 0; icor < numcor; icor++) {
            if (ACFcount[icor] > 0) {
                fprintf(output_file, " %f  \t %f\n", icor * num_dtau * DT, cor[icor] / (double)ACFcount[icor]);
            }
        }
        fclose(output_file);
    }
}

/**
 * @brief Calculates and outputs velocity profiles.
 * @param reset 0: initialize, 1: update, 2: finalize and write output.
 * @param time Current simulation time.
 */
void velocityprofile(int reset, double time)
{
    FILE* output_file;
    char filename[256]; // Buffer for filename

    // Construct filename based on time
    // Original path was hardcoded, now relative to data/
    snprintf(filename, sizeof(filename), "data/velocity_profile=%.0f_%.0f.dat", floor(time), (time - floor(time)) * 1000);

    if (reset == 0) { // Initialize
        ybin = 100; // Number of bins in y-direction
        vel_sum = (double*)calloc(ybin, sizeof(double));
        counter = (int*)calloc(ybin, sizeof(int));
        // File is opened in reset==2, not here.
    } else if (reset == 1) { // Update
        dbin = box_y / (double)ybin;
        for (int i = 0; i < Npart; i++) {
            // Only consider particles within a certain z-range (e.g., center of the box)
            if (Rz[i] - (double)L_z / 2.0 <= 1.0 && Rz[i] - (double)L_z / 2.0 >= -1.0) {
                int ibin = (int)(Ry[i] / dbin);
                if (ibin >= 0 && ibin < ybin) {
                    vel_sum[ibin] += Vz[i];
                    counter[ibin] += 1;
                }
            }
        }
    } else { // Finalize and write output
        output_file = fopen(filename, "w");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        for (int ibin = 0; ibin < ybin; ibin++) {
            if (counter[ibin] > 0) {
                fprintf(output_file, "%f \t %f \n", ibin * dbin, vel_sum[ibin] / (double)counter[ibin]);
            }
            else {
                fprintf(output_file, "%f \t %f \n", ibin * dbin, 0.0); // Output 0 if no particles in bin
            }
        }
        fclose(output_file);
        free_1D_array(vel_sum);
        free_1D_array(counter);
    }
}

/**
 * @brief Tracks particles in a slit geometry and outputs their y-distribution.
 * @param reset 0: initialize, 1: update and write output.
 * @param time Current simulation time.
 */
void slit_tracker(int reset, double time)
{
    FILE* output_file;
    char filename[256]; // Buffer for filename

    // Construct filename based on time
    // Original path was hardcoded, now relative to data/
    snprintf(filename, sizeof(filename), "data/tracker_time_profile_1=%.0f_%.0f.dat", floor(time), (time - floor(time)) * 1000);

    if (reset == 0) { // Initialize
        ybin = 100; // Number of bins in y-direction
        y_sum = (double*)calloc(ybin, sizeof(double));
    } else if (reset == 1) { // Update and write output
        dbin = (double)L_y / (double)ybin;
        for (int i = 0; i < num_track; i++) {
            int ibin = (int)(Ry[tracker[i]] / dbin);
            if (ibin >= 0 && ibin < ybin) { // Ensure bin index is valid
                y_sum[ibin] += 1.0;
            }
        }
        output_file = fopen(filename, "a");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        for (int ibin = 0; ibin < ybin; ibin++) {
            fprintf(output_file, "%f \t %f \n", ibin * dbin, y_sum[ibin]);
            y_sum[ibin] = 0.0; // Reset for next interval
        }
        fclose(output_file);
    }
}

/**
 * @brief Creates and outputs concentration profiles for different regions.
 * This function is used for analyzing particle distribution in specific areas.
 * @param reset 0: initialize, 1: update, 2: finalize and write output.
 * @param time Current simulation time.
 */
void profile_maker(int reset, double time)
{
    FILE* output_file_center;
    FILE* output_file_edge;
    FILE* output_file_total;
    char filename_center[256], filename_edge[256], filename_total[256];

    // Construct filenames based on time
    // Original path was hardcoded, now relative to data/
    snprintf(filename_center, sizeof(filename_center), "data/tracker_time_center=%.0f_%.0f.dat", floor(time), (time - floor(time)) * 1000);
    snprintf(filename_edge, sizeof(filename_edge), "data/tracker_time_edge=%.0f_%.0f.dat", floor(time), (time - floor(time)) * 1000);
    snprintf(filename_total, sizeof(filename_total), "data/tracker_time_total=%.0f_%.0f.dat", floor(time), (time - floor(time)) * 1000);

    if (reset == 0) { // Initialize
        zbin = 200; // Number of bins in z-direction
        z_sum = (double*)calloc(zbin, sizeof(double));
        z_sum1 = (double*)calloc(zbin, sizeof(double));
    } else if (reset == 1) { // Update
        dbin = box_z / (double)zbin;
        for (int i = 0; i < Npart; i++) {
            if (!(Ry[i] < beta * box_y || Ry[i] >= (1.0 - beta) * box_y)) { // Center region
                int ibin = (int)(Rz[i] / dbin);
                if (ibin >= 0 && ibin < zbin) {
                    z_sum[ibin] += 1.0;
                }
            } else { // Edge region
                int ibin2 = (int)(Rz[i] / dbin);
                if (ibin2 >= 0 && ibin2 < zbin) {
                    z_sum1[ibin2] += 1.0;
                }
            }
        }

        output_file_center = fopen(filename_center, "a");
        output_file_edge = fopen(filename_edge, "a");
        output_file_total = fopen(filename_total, "a");

        if (output_file_center == NULL || output_file_edge == NULL || output_file_total == NULL) {
            fprintf(stderr, "Error: Could not open profile maker files for writing.\n", filename);
            return;
        }

        for (int ibin = 0; ibin < zbin; ibin++) {
            // Normalize by volume of the region
            double volume_center = box_z * (1 - 2 * beta) * box_y * box_x;
            double volume_edge = box_z * (2 * beta) * box_y * box_x;
            double volume_total = box_z * box_y * box_x;

            fprintf(output_file_center, "%f \t %f \n", ibin * dbin, (volume_center > 0) ? (double)zbin * z_sum[ibin] / volume_center : 0.0);
            fprintf(output_file_edge, "%f \t %f \n", ibin * dbin, (volume_edge > 0) ? (double)zbin * z_sum1[ibin] / volume_edge : 0.0);
            fprintf(output_file_total, "%f \t %f \n", ibin * dbin, (volume_total > 0) ? (double)zbin * (z_sum1[ibin] + z_sum[ibin]) / volume_total : 0.0);

            z_sum[ibin] = 0.0; // Reset for next interval
            z_sum1[ibin] = 0.0;
        }
        fclose(output_file_center);
        fclose(output_file_edge);
        fclose(output_file_total);
    } else { // Finalize (no specific action needed for reset==2, memory freed in main)
        free_1D_array(z_sum);
        free_1D_array(z_sum1);
    }
}

/**
 * @brief Calculates and outputs unary velocity autocorrelation function (ACF).
 * @param reset 1: initialize, 0: update, 2: finalize and write output.
 * @param DT Time step.
 */
void velocity_correlator(int reset, double DT)
{
    FILE* output_file;
    char filename[256]; // Buffer for filename

    // Construct filename based on DT
    snprintf(filename, sizeof(filename), "data/vel_Corr_dt=%.3f.dat", DT);

    if (reset == 1) { // Initialize
        frame = 0;
        output_file = fopen(filename, "w");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        fprintf(output_file, "velocity\t\t time \n");
        fclose(output_file);
    } else if (reset == 0) { // Update
        int maxcor_val = numcor - 1;
        if (frame < numcor) maxcor_val = frame;
        int curframe = frame % numcor;

        for (int i = 0; i < Npart; i++) {
            storex[i][curframe] = Vx[i];
            storey[i][curframe] = Vy[i];
            storez[i][curframe] = Vz[i];

            for (int icor = 0; icor <= maxcor_val; icor++) {
                int icorframe = (curframe - icor + numcor) % numcor;
                cor[icor] += Vx[i] * storex[i][icorframe] + Vy[i] * storey[i][icorframe] + Vz[i] * storez[i][icorframe];
                ACFcount[icor] += 1;
            }
        }
        frame++;
    } else { // Finalize and write output
        double final_diffusion_coeff = 0.0;
        output_file = fopen(filename, "a");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        for (int icor = 0; icor < numcor; icor++) {
            if (ACFcount[icor] > 0) {
                double current_acf = cor[icor] / (double)ACFcount[icor];
                fprintf(output_file, " %f  \t %f\n", icor * num_dtau * DT, current_acf);
                if (icor < numcor - 1) {
                    double next_acf = cor[icor + 1] / (double)ACFcount[icor + 1];
                    final_diffusion_coeff += (current_acf + next_acf) * DT / 2.0; // Trapezoidal rule for integration
                }
            }
        }
        fprintf(output_file, "Diffusion Coefficient = %f\n", final_diffusion_coeff / 3.0);
        fclose(output_file);
        printf("\n***Diffusion Coefficient = %f***\n", final_diffusion_coeff / 3.0);
    }
}

/**
 * @brief Calculates and outputs dual velocity autocorrelation function (ACF) for two species.
 * This is used for mutual diffusion coefficients.
 * @param reset 1: initialize, 0: update, 2: finalize and write output.
 * @param massA Mass of foreign particles.
 */
void dual_velocity_correlator(int reset, double massA)
{
    FILE* output_file;
    char filename[256]; // Buffer for filename

    // Construct filename based on concentration and massA
    snprintf(filename, sizeof(filename), "data/(N)_Diffusion_concA=%.3f_m=%.1f.dat", conc, massA);

    if (reset == 1) { // Initialize
        frame = 0;
        output_file = fopen(filename, "w");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        fprintf(output_file, "velocity\t\t time \n");
        fclose(output_file);
    } else if (reset == 0) { // Update
        int maxcor_val = numcor - 1;
        if (frame < numcor) maxcor_val = frame;
        int curframe = frame % numcor;

        // For SRD particles
        for (int i = 0; i < Nsrd; i++) {
            storex[i][curframe] = Vx[i];
            storey[i][curframe] = Vy[i];
            storez[i][curframe] = Vz[i];
            for (int icor = 0; icor <= maxcor_val; icor++) {
                int icorframe = (curframe - icor + numcor) % numcor;
                cor[icor] += Vx[i] * storex[i][icorframe] + Vy[i] * storey[i][icorframe] + Vz[i] * storez[i][icorframe];
                ACFcount[icor] += 1;
            }
        }
        // For foreign particles
        for (int i = Nsrd; i < Npart; i++) {
            storex1[i][curframe] = Vx[i];
            storey1[i][curframe] = Vy[i];
            storez1[i][curframe] = Vz[i];
            for (int icor = 0; icor <= maxcor_val; icor++) {
                int icorframe = (curframe - icor + numcor) % numcor;
                cor1[icor] += Vx[i] * storex1[i][icorframe] + Vy[i] * storey1[i][icorframe] + Vz[i] * storez1[i][icorframe];
                ACFcount1[icor] += 1;
            }
        }
        frame++;
    } else { // Finalize and write output
        double helper_srd = 0.0;
        double helper_foreign = 0.0;
        output_file = fopen(filename, "a");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        for (int icor = 0; icor < numcor; icor++) {
            double current_acf_srd = (ACFcount[icor] > 0) ? cor[icor] / (double)ACFcount[icor] : 0.0;
            double current_acf_foreign = (ACFcount1[icor] > 0) ? cor1[icor] / (double)ACFcount1[icor] : 0.0;
            fprintf(output_file, " %f  \t  %f\n", icor * num_dtau * DT, current_acf_srd, current_acf_foreign);

            if (icor < numcor - 1) {
                double next_acf_srd = (ACFcount[icor + 1] > 0) ? cor[icor + 1] / (double)ACFcount[icor + 1] : 0.0;
                double next_acf_foreign = (ACFcount1[icor + 1] > 0) ? cor1[icor + 1] / (double)ACFcount1[icor + 1] : 0.0;
                helper_srd += (current_acf_srd + next_acf_srd) * DT / 2.0;
                helper_foreign += (current_acf_foreign + next_acf_foreign) * DT / 2.0;
            }
        }
        fprintf(output_file, "Diffusion Coefficient : %f \t %f\n", helper_srd / 3.0, helper_foreign / 3.0);
        fclose(output_file);
        printf("\nDiff SRD: %f \t \t ::: A: %f\n", helper_srd / 3.0, helper_foreign / 3.0);
    }
}

/**
 * @brief Calculates and outputs terminal velocity for diffusion coefficient calculation.
 * This method is based on applying an external force and measuring terminal velocity.
 * @param reset 1: initialize, 0: update, 2: finalize and write output.
 * @param t Current simulation time.
 */
void new_diff(int reset, double t)
{
    FILE* output_file;
    char filename[256]; // Buffer for filename

    // Construct filename based on external force 'g'
    snprintf(filename, sizeof(filename), "data/g=%.3f.dat", g);

    if (reset == 1) { // Initialize
        average = 0.0;
        cc = 0;
        frame = 0;
        average1 = 0.0;
        cc1 = 0;
        output_file = fopen(filename, "w");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        fprintf(output_file, "velocity\t\t time \n");
        fclose(output_file);
    } else if (reset == 0) { // Update
        double U_srd = 0.0;
        double U_foreign = 0.0;
        for (int i = 0; i < Nsrd; i++) {
            U_srd += Vz[i];
        }
        for (int i = Nsrd; i < Npart; i++) {
            U_foreign += Vz[i];
        }
        average = (average * cc + U_srd) / (++cc);
        average1 = (average1 * cc1 + U_foreign) / (++cc1);

        output_file = fopen(filename, "a");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        fprintf(output_file, "%f \t %f    \n", U_srd, t); // Output SRD particle velocity sum
        fclose(output_file);
    } else { // Finalize and print average terminal velocity
        printf("***Average SRD Terminal Velocity: %f****\n", average / (double)Nsrd);
        printf("***Average Foreign Terminal Velocity: %f****\n", average1 / (double)NA);
    }
}

/**
 * @brief Calculates and outputs Green-Kubo relation for mutual diffusion coefficient.
 * This function is named 'something' in the original code.
 * @param reset 1: initialize, 0: update, 2: finalize and write output.
 * @param massA Mass of foreign particles.
 */
void something(int reset, double massA)
{
    double ratio = 0.0;
    if (conc != 0) ratio = (double)Nsrd / (massA * (double)NA);

    FILE* output_file;
    char filename[256]; // Buffer for filename

    // Construct filename based on concentration and massA
    snprintf(filename, sizeof(filename), "data/(N)_Diffusion_concA=%.3f_m=%.1f.dat", conc, massA);

    if (reset == 1) { // Initialize
        frame = 0;
        output_file = fopen(filename, "w");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        fprintf(output_file, "velocity\t\t time \n");
        fclose(output_file);
    } else if (reset == 0) { // Update
        double sum_vx_srd = 0.0, sum_vy_srd = 0.0, sum_vz_srd = 0.0;
        double sum_vx_foreign = 0.0, sum_vy_foreign = 0.0, sum_vz_foreign = 0.0;

        int maxcor_val = numcor - 1;
        if (frame < numcor) maxcor_val = frame;
        int curframe = frame % numcor;

        // Calculate sum of velocities for SRD particles
        for (int i = 0; i < Nsrd; i++) {
            sum_vx_srd += Vx[i];
            sum_vy_srd += Vy[i];
            sum_vz_srd += Vz[i];
            storex1[i][curframe] = Vx[i];
            storey1[i][curframe] = Vy[i];
            storez1[i][curframe] = Vz[i];
            for (int icor = 0; icor <= maxcor_val; icor++) {
                int icorframe = (curframe - icor + numcor) % numcor;
                cor1[icor] += Vx[i] * storex1[i][icorframe] + Vy[i] * storey1[i][icorframe] + Vz[i] * storez1[i][icorframe];
                ACFcount1[icor] += 1;
            }
        }
        // Calculate sum of velocities for foreign particles
        for (int i = Nsrd; i < Npart; i++) {
            sum_vx_foreign += Vx[i];
            sum_vy_foreign += Vy[i];
            sum_vz_foreign += Vz[i];
            storex2[i][curframe] = Vx[i];
            storey2[i][curframe] = Vy[i];
            storez2[i][curframe] = Vz[i];
            for (int icor = 0; icor <= maxcor_val; icor++) {
                int icorframe = (curframe - icor + numcor) % numcor;
                cor2[icor] += Vx[i] * storex2[i][icorframe] + Vy[i] * storey2[i][icorframe] + Vz[i] * storez2[i][icorframe];
                ACFcount2[icor] += 1;
            }
        }

        // Calculate relative velocity difference (SRD - Foreign)
        if (Nsrd > 0 && NA > 0) {
            storex[0][curframe] = sum_vx_srd / (double)Nsrd - sum_vx_foreign / (double)NA;
            storey[0][curframe] = sum_vy_srd / (double)Nsrd - sum_vy_foreign / (double)NA;
            storez[0][curframe] = sum_vz_srd / (double)Nsrd - sum_vz_foreign / (double)NA;
        } else {
            storex[0][curframe] = 0.0;
            storey[0][curframe] = 0.0;
            storez[0][curframe] = 0.0;
        }


        for (int icor = 0; icor <= maxcor_val; icor++) {
            int icorframe = (curframe - icor + numcor) % numcor;
            // Correlate SRD sum velocity with the relative velocity difference
            cor[icor] += (sum_vx_srd / (double)Nsrd) * storex[0][icorframe] +
                         (sum_vy_srd / (double)Nsrd) * storey[0][icorframe] +
                         (sum_vz_srd / (double)Nsrd) * storez[0][icorframe];
            ACFcount[icor] += 1;
        }
        frame++;
    } else { // Finalize and write output
        double helper_mixture = 0.0;
        double helper_srd = 0.0;
        double helper_foreign = 0.0;

        output_file = fopen(filename, "a");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        for (int icor = 0; icor < numcor; icor++) {
            double current_cor_mixture = (ACFcount[icor] > 0) ? cor[icor] / (double)ACFcount[icor] : 0.0;
            double current_cor_srd = (ACFcount1[icor] > 0) ? cor1[icor] / (double)ACFcount1[icor] : 0.0;
            double current_cor_foreign = (ACFcount2[icor] > 0) ? cor2[icor] / (double)ACFcount2[icor] : 0.0;

            fprintf(output_file, " %f  \t  %f  \n", icor * num_dtau * DT, current_cor_mixture);

            // Simpson's rule for integration (original code used this pattern)
            if (icor == 0 || icor == numcor - 1) { // First and last point
                helper_mixture += DT / 3.0 * current_cor_mixture;
                helper_srd += DT / 3.0 * current_cor_srd;
                helper_foreign += DT / 3.0 * current_cor_foreign;
            } else if (icor % 2 == 0) { // Even points
                helper_mixture += 2.0 * DT / 3.0 * current_cor_mixture;
                helper_srd += 2.0 * DT / 3.0 * current_cor_srd;
                helper_foreign += 2.0 * DT / 3.0 * current_cor_foreign;
            } else { // Odd points
                helper_mixture += 4.0 * DT / 3.0 * current_cor_mixture;
                helper_srd += 4.0 * DT / 3.0 * current_cor_srd;
                helper_foreign += 4.0 * DT / 3.0 * current_cor_foreign;
            }
        }
        // Calculate mutual diffusion coefficient
        double mutual_diffusion_coeff = (1.0 + ratio) * (1.0 + ratio) * helper_mixture / ((double)Npart * ratio * massA * 3.0);
        fprintf(output_file, "Diffusion Coefficient : %f \n", mutual_diffusion_coeff);
        fclose(output_file);
        printf("\n**Diff mixture: %f*** ; Diffusion SRD : %f*** ; Diffusion foreign : %f***\n", mutual_diffusion_coeff, helper_srd / 3.0, helper_foreign / 3.0);
    }
}
