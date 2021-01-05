#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> // For snprintf

#include "constants.h"
#include "initialise.h"
#include "functions.h"
#include "variables.h"
#include "memory.h" // For memory allocation functions

/**
 * @brief Manages homogeneous particle distribution and counts particles
 * within a specific region for axial dispersion method (CONV_DIFF == 2).
 * @param reset 0: initialize, 1: update, 2: finalize and write output.
 * @param tx Current time modulo t_interval.
 */
void homogenous(int reset, double tx)
{
    if (reset == 0) { // Initialize
        arrays = (int*)calloc((int)(t_interval / dt), sizeof(int));
    } else if (reset == 1) { // Update
        for (int i = Nsrd; i < Npart; i++) { // Iterate over foreign particles
            // Check if particle is within a specific z-range (e.g., measurement region)
            if (Rz_artificial[i] < (1 + alpha) / 2.0 * box_z + 0.5 && Rz_artificial[i] > (1 + alpha) / 2.0 * box_z - 0.5) {
                if ((int)tx >= 0 && (int)tx < (int)(t_interval / dt)) { // Ensure tx is a valid index
                    arrays[(int)tx] += 1;
                }
            }
        }
    } else { // Finalize and write output
        FILE* output_file;
        char filename[256];
        snprintf(filename, sizeof(filename), "data/Homogenous_3.dat"); // Hardcoded filename

        output_file = fopen(filename, "w");
        if (output_file == NULL) {
            fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
            return;
        }
        for (int i = 0; i < (int)(t_interval / dt); i++) {
            fprintf(output_file, " %f  \t %d\n", (double)(i * dt), arrays[i]);
        }
        fclose(output_file);
        free_1D_array(arrays);
    }
}

/**
 * @brief Updates z-velocity bins for flow profile generation.
 * This function is used in conjunction with velocity_changer to impose a flow profile.
 * @param reset 0: initialize, 1: update.
 * @param total_time Current simulation time.
 */
void update_velz_bin(int reset, int total_time)
{
    if (reset == 0) { // Initialize
        velz_bin = L_y * 5; // Number of bins in y-direction for velocity profile
        velz_meas = (double*)calloc(velz_bin, sizeof(double));
        velz_flow = (double*)calloc(velz_bin, sizeof(double));
        N_average = (double*)calloc(velz_bin, sizeof(double));

        // Initialize velz_flow with Gaussian random numbers
        for (int i = 0; i < velz_bin; i++) {
            velz_flow[i] = sqrt(kTo) * neargauss();
        }
        N_recur = 0.0; // Reset N_recur (global variable, purpose unclear without more context)
    } else if (reset == 1) { // Update
        double dbin = box_y / (double)velz_bin;
        double N_temp = 0.0; // Temporary particle count in measurement region
        double N_recur_local = 0.0; // Local N_recur

        for (int i = 0; i < Npart; i++) {
            // Check if particle is within a specific z-range for measurement
            if (Rz[i] < box_z * (1.0 + alpha) / 2.0 + 0.5 && Rz[i] > box_z * (1.0 + alpha) / 2.0 - 0.5) {
                N_temp += 1.0;
                N_recur_local += 1.0;
                int y_bin = (int)(Ry[i] / dbin);
                if (y_bin >= 0 && y_bin < velz_bin) {
                    velz_meas[y_bin] = (velz_meas[y_bin] * N_average[y_bin] + Vz[i]) / (N_average[y_bin] + 1.0);
                    N_average[y_bin] += 1.0;
                }
            }
        }

        // Update velz_flow based on measured velocities and an exponential decay
        for (int i = 0; i < velz_bin; i++) {
            if (N_recur_local > 0) { // Avoid division by zero
                velz_flow[i] = velz_flow[i] * exp(-sqrt(10 * N_temp / N_recur_local)) + velz_meas[i] * (1 - exp(-sqrt(10 * N_temp / N_recur_local)));
            } else {
                velz_flow[i] = velz_flow[i]; // No update if no particles in region
            }
            velz_meas[i] = 0.0; // Reset measured velocities
            N_average[i] = 0.0; // Reset average counts
        }
    }
}

/**
 * @brief Changes the velocity of a particle based on the imposed flow profile.
 * This function is called for individual particles.
 * @param i Index of the particle.
 */
void velocity_changer(int i)
{
    int cell_y_bin = (int)(Ry[i] * (double)velz_bin / box_y);
    if (cell_y_bin < 0) cell_y_bin = 0;
    if (cell_y_bin >= velz_bin) cell_y_bin = velz_bin - 1;

    double Ry_temp_high = (cell_y_bin == velz_bin - 1) ? box_y : box_y / (double)velz_bin * (cell_y_bin + 1);
    double Ry_temp_low = box_y / (double)velz_bin * cell_y_bin;

    // Linear interpolation for Vz
    if (cell_y_bin == velz_bin - 1) { // Last bin, use value from first bin for interpolation (periodic-like)
        Vz[i] = velz_flow[0] - (Ry_temp_high - Ry[i]) * (velz_flow[0] - velz_flow[cell_y_bin]) / (Ry_temp_high - Ry_temp_low);
    } else {
        Vz[i] = velz_flow[cell_y_bin + 1] - (Ry_temp_high - Ry[i]) * (velz_flow[cell_y_bin + 1] - velz_flow[cell_y_bin]) / (Ry_temp_high - Ry_temp_low);
    }

    // Assign new random velocities for Vx and Vy (thermal motion)
    Vx[i] = sqrt(kTo) * neargauss();
    Vy[i] = sqrt(kTo) * neargauss();
}
