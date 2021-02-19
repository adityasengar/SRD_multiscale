#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "constants.h"
#include "initialise.h"
#include "functions.h"
#include "variables.h"
#include "memory.h" // For memory allocation functions

// Global variable for cell-linked-list (c is already a global in original, but not used consistently)
// int c = 0; // Removed, as it's not used consistently and can be local if needed.

/**
 * @brief Calculates the modulo of a double r with an integer a.
 * Handles both positive and negative r correctly.
 * @param r The double number.
 * @param a The integer divisor.
 * @return The result of r modulo a.
 */
double mod(double r, int a)
{
    double x;
    if (r > 0.0) x = r - (double)(a * (int)(r / a));
    else x = r - (double)((a * (int)((r / a) - 1)));
    return x;
}

/**
 * @brief Creates a cell-linked-list for efficient particle neighbor searching.
 * Includes logic for random grid shift and particle combining/splitting
 * for axial dispersion methods.
 *
 * head[cellx][celly][cellz] is the first particle in the cell.
 * list[i] is the next particle in the same cell.
 * If list[i] == -1, the last particle in the cell is reached.
 */
void createlinkedlist()
{
    int cellx, celly, cellz, i;
    int count_pair = 0; // Counter for combined particles

    // --- Particle Combining/Splitting Logic for Axial Dispersion (CONV_DIFF_mass == 1) ---
    // This section handles the creation/destruction of particles based on their 'activate' status
    // and position, specifically for the axial dispersion method with different mass particles.
    if (CONV_DIFF_mass == 1 && walls == 1 && CONV_DIFF == 2 && t >= tmin && (double)((int)(t / DT) % (int)(t_interval / DT)) < (t_step / DT))
    {
        // Count particles marked for combining
        int marked_for_combining = 0;
        for (cellx = 0; cellx < L_x; cellx++) {
            for (celly = 0; celly < L_y + 1; celly++) {
                for (cellz = 0; cellz < (int)(alpha * L_z); cellz++) {
                    int j = head_mass[cellx][celly][cellz];
                    while (j >= 0) {
                        if (activate[j][0] == -1) marked_for_combining++;
                        j = list_mass[j];
                    }
                }
            }
        }
        printf("\nNumber of particles that are marked for combining: %d\n", marked_for_combining);

        // Perform combining of particles
        for (cellx = 0; cellx < L_x; cellx++) {
            for (celly = 0; celly < L_y + 1; celly++) {
                for (cellz = 0; cellz < (int)(alpha * L_z); cellz++) {
                    int j = head_mass[cellx][celly][cellz];
                    while (j >= 0) {
                        if (list_mass[j] >= 0 && activate[j][1] == 0) {
                            // Combine two particles into a new one
                            Rx[Npart] = (Rx[j] + Rx[list_mass[j]]) / 2.0;
                            Ry[Npart] = (Ry[j] + Ry[list_mass[j]]) / 2.0;
                            Rz[Npart] = (Rz[j] + Rz[list_mass[j]]) / 2.0;
                            Rz_artificial[Npart] = Rz[Npart];
                            Vx[Npart] = (Vx[j] + Vx[list_mass[j]]) / 2.0;
                            Vy[Npart] = (Vy[j] + Vy[list_mass[j]]) / 2.0;
                            Vz[Npart] = (Vz[j] + Vz[list_mass[j]]) / 2.0;

                            activate[list_mass[j]][1] = -1; // Mark original particles for removal
                            activate[j][1] = -1;

                            activate[Npart][0] = -1; // Mark new particle as combined
                            activate[Npart][1] = 0;

                            Npart++;
                            count_pair++;
                        } else {
                            activate[j][1] = 0;
                            activate[j][0] = 0;
                        }
                        // Move to the next pair or particle
                        if (list_mass[j] >= 0 && list_mass[list_mass[j]] >= 0)
                            j = list_mass[list_mass[j]];
                        else
                            j = -1;
                    }
                }
            }
        }
        printf("!!!%d!!!", count_pair);
        NA = count_pair; // Update number of foreign particles (A)

        // Remove particles marked for removal by shifting array elements
        i = 0;
        while (i < Npart) {
            if (activate[i][1] == -1 && activate[i][0] == -1) {
                for (int k = i; k < Npart - 1; k++) {
                    Rx[k] = Rx[k + 1];
                    Ry[k] = Ry[k + 1];
                    Rz[k] = Rz[k + 1];
                    Rz_artificial[k] = Rz_artificial[k + 1];
                    Vx[k] = Vx[k + 1];
                    Vy[k] = Vy[k + 1];
                    Vz[k] = Vz[k + 1];
                    activate[k][0] = activate[k + 1][0];
                    activate[k][1] = activate[k + 1][1];
                }
                Nsrd--;
                Npart--;
            } else {
                i++;
            }
        }
    }
    // --- End of Particle Combining/Splitting Logic ---

    // Reset head array for linked list
    for (cellx = 0; cellx < L_x; cellx++) {
        for (celly = 0; celly < L_y + 1; celly++) {
            for (cellz = 0; cellz < L_z + 1; cellz++) { // L_z + 1 for open boundary
                head[cellx][celly][cellz] = -1;
            }
        }
    }

    // Reset list array
    for (i = 0; i < Npart; i++) {
        list[i] = -1;
    }

    // Apply random grid shift and populate linked list
    if (doshift == 1) {
        shiftx = (double)rand() / (double)RAND_MAX;
        shifty = (double)rand() / (double)RAND_MAX;
        shiftz = (double)rand() / (double)RAND_MAX;
        for (i = 0; i < Npart; i++) {
            cellx = (int)(Rx[i] + shiftx);
            celly = (int)(Ry[i] + shifty);
            cellz = (int)(Rz[i] + shiftz);

            cellx = cellx % L_x;
            if (open == 0) { // Periodic boundaries in z-direction
                cellz = cellz % L_z;
            }
            if (walls == 0) { // Periodic boundaries in y
                celly = celly % L_y;
            }
            // For solid walls, celly can go beyond L_y, handled by nearwalls and ghost particles

            list[i] = head[cellx][celly][cellz];
            head[cellx][celly][cellz] = i;
        }
    }
    else { // No grid-shift applied
        shiftx = 0.0;
        shifty = 0.0;
        shiftz = 0.0;
        for (i = 0; i < Npart; i++) {
            cellx = (int)(Rx[i]);
            celly = (int)(Ry[i]);
            cellz = (int)(Rz[i]);

            list[i] = head[cellx][celly][cellz];
            head[cellx][celly][cellz] = i;
        }
    }

    // Particle splitting for CONV_DIFF_mass == 1 (if conditions met)
    if (CONV_DIFF_mass == 1 && walls == 1 && CONV_DIFF == 2 && t >= tmin && (double)((int)(t / DT) % (int)(t_interval / DT)) >= 2900.0 && (double)((int)(t / DT) % (int)(t_interval / DT)) < 2901.0)
    {
        int Npart_temp = Npart;
        int Nsrd_temp = Nsrd;
        for (i = 0; i < Npart_temp; i++) {
            if (i >= Nsrd_temp) { // Only split foreign particles (A)
                Rx[Npart] = Rx[i];
                Ry[Npart] = Ry[i];
                Rz[Npart] = Rz[i];

                // Split velocities with Gaussian noise
                double Vx1 = neargauss();
                double Vx2 = neargauss();
                double Vy1 = neargauss();
                double Vy2 = neargauss();
                double Vz1 = neargauss();
                double Vz2 = neargauss();

                Vx[i] = Vx[i] + Vx1 - (Vx1 + Vx2) / 2.0;
                Vx[Npart] = Vx[i] + Vx2 - (Vx1 + Vx2) / 2.0;
                Vy[i] = Vy[i] + Vy1 - (Vy1 + Vy2) / 2.0;
                Vy[Npart] = Vy[i] + Vy2 - (Vy1 + Vy2) / 2.0;
                Vz[i] = Vz[i] + Vz1 - (Vz1 + Vz2) / 2.0;
                Vz[Npart] = Vz[i] + Vz2 - (Vz1 + Vz2) / 2.0;

                activate[i][0] = 0; // Mark as active
                activate[i][1] = 0;
                activate[Npart][0] = 0;
                activate[Npart][1] = 0;
                Npart++;
                NA--;
                Nsrd += 2; // Two new SRD particles are created from one foreign particle
            } else {
                activate[i][0] = 0;
                activate[i][1] = 0;
            }
        }
        printf("Coalescence : Nsrd, NA, Ntot :%f, %d, %d, %d;;;\n\n", t, Nsrd, NA, Npart);

        // Reallocate and reinitialize mass-related linked list arrays after particle count change
        activate = int_2D_array(2 * Npart, 4); // Assuming activate needs to be resized
        list_mass = int_1D_array(Npart);
        head_mass = int_3D_array(LBOX_MAX, LBOX_MAX + 1, LBOX_MAX); // Use LBOX_MAX from constants.h

        for (cellx = 0; cellx < L_x; cellx++) {
            for (celly = 0; celly < L_y + 1; celly++) {
                for (cellz = 0; cellz < L_z; cellz++) {
                    head_mass[cellx][celly][cellz] = -1;
                }
            }
        }
        for (i = 0; i < Npart; i++) {
            list_mass[i] = -1;
        }
    }
}

/**
 * @brief Updates velocities of SRD particles due to external fields (gravity).
 */
void updatevelocities()
{
    for (int i = 0; i < Nsrd; i++) { // SRD particles
        Vz[i] += g * DT;
    }
    for (int i = Nsrd; i < Npart; i++) { // Foreign particles
        Vz[i] += g * DT; // Assuming foreign particles also experience gravity
    }
}

/**
 * @brief Returns the minimum of three double numbers.
 * @param r1 First number.
 * @param r2 Second number.
 * @param r3 Third number.
 * @return The smallest of the three numbers.
 */
double min3(double r1, double r2, double r3)
{
    double lowest = r1;
    if (r2 < lowest) lowest = r2;
    if (r3 < lowest) lowest = r3;
    return lowest;
}

/**
 * @brief Updates particle positions and handles boundary conditions.
 * Includes logic for solid walls (bounce-back) and catalytic surface interactions.
 */
void updatepositions()
{
    int i;
    double xprev, yprev, zprev, xnew, ynew, znew;
    double time_to_wall_y; // Renamed from 'timey' to avoid conflict with global 'timey'

    // Reallocate 'activate' array if needed (original code had this logic)
    // This part of the original code was commented out or had conditional allocation
    // if((double)((int)t%(int)t_interval)==t_interval-1.0) {
    //     activate=int_2D_array(Npart*2,4);
    // }

    for (i = 0; i < Npart; i++) { // Update SRD positions
        xprev = Rx[i];
        yprev = Ry[i];
        zprev = Rz[i];
        xnew = xprev + Vx[i] * DT;
        ynew = yprev + Vy[i] * DT;
        znew = zprev + Vz[i] * DT;
        Rz_artificial[i] += Vz[i] * DT; // Update artificial z-position

        // Handle boundary conditions
        if (xnew < 0.0 || xnew >= box_x || ynew < 0.0 || ynew >= box_y || znew < 0.0 || znew >= box_z) {
            // Periodic boundaries for x
            if (xnew >= box_x || xnew < 0) {
                xnew = mod(xnew, L_x);
            }
            // Periodic or open boundaries for z
            if (open == 0) { // Periodic boundaries in z
                if (znew < 0 || znew >= box_z) {
                    znew = mod(znew, L_z);
                }
            } else { // Open system in z, particles leaving are removed
                if (znew < 0) znew = 0.0;
                if (znew >= box_z) znew = box_z - 0.0001; // Slightly less than box_z
            }


            if (walls == 0) { // Periodic boundaries for y
                if (ynew < 0.0 || ynew >= box_y) {
                    ynew = mod(ynew, L_y);
                }
            } else { // Solid walls (bounce-back) for y
                if (ynew < 0.0) { // Collision with bottom wall (y=0)
                    time_to_wall_y = ynew / (Vy[i]); // Time to reach the wall
                    double xnew_surface = xnew - time_to_wall_y * Vx[i];
                    double ynew_surface = ynew - time_to_wall_y * Vy[i];
                    double znew_surface = znew - time_to_wall_y * Vz[i];

                    // Ensure surface coordinates are within periodic boundaries for x and z
                    if (xnew_surface >= box_x || xnew_surface < 0) xnew_surface = mod(xnew_surface, L_x);
                    if (open == 0) { // Periodic boundaries in z
                        if (znew_surface < 0 || znew_surface >= box_z) znew_surface = mod(znew_surface, L_z);
                    } else { // Open system in z
                        if (znew_surface < 0) znew_surface = 0.0;
                        if (znew_surface >= box_z) znew_surface = box_z - 0.0001;
                    }


                    // Catalytic surface interaction logic
                    // This region (z: 4.0-6.0) is hardcoded as catalytic
                    if (znew_surface > 4.0 && znew_surface < 6.0) {
                        double random_val = (double)rand() / (double)RAND_MAX;
                        if (i < Nsrd) { // SRD particle (potential reactant A)
                            if (random_val < p_adsA * theta) { // Adsorption of A
                                catalyst[i][0] = 0; // Mark as SRD particle
                                catalyst[i][1] = 1; // Mark as adsorbed A
                                Rx[i] = xnew_surface;
                                Ry[i] = ynew_surface;
                                Rz[i] = znew_surface;
                                continue; // Skip bounce-back for adsorbed particles
                            }
                        } else { // Foreign particle (potential reactant B or product)
                            if (random_val < p_adsB * theta) { // Adsorption of B
                                catalyst[i][0] = 3; // Mark as foreign particle
                                catalyst[i][1] = 2; // Mark as adsorbed B
                                Rx[i] = xnew_surface;
                                Ry[i] = ynew_surface;
                                Rz[i] = znew_surface;
                                continue; // Skip bounce-back for adsorbed particles
                            }
                        }
                    }

                    // If not adsorbed, perform bounce-back
                    if (catalyst[i][1] == 0 || catalyst[i][1] == 3) { // Only bounce back if not adsorbed
                        xnew = xnew - 2.0 * time_to_wall_y * Vx[i];
                        ynew = ynew - 2.0 * time_to_wall_y * Vy[i];
                        znew = znew - 2.0 * time_to_wall_y * Vz[i];
                        Rz_artificial[i] -= 2.0 * time_to_wall_y * Vz[i];
                    }
                    Vx[i] = -Vx[i]; // Invert velocity
                    Vy[i] = -Vy[i];
                    Vz[i] = -Vz[i];
                }
                if (ynew >= box_y) { // Collision with top wall (y=box_y)
                    time_to_wall_y = (ynew - box_y) / (Vy[i]);
                    xnew = xnew - 2.0 * time_to_wall_y * Vx[i];
                    ynew = ynew - 2.0 * time_to_wall_y * Vy[i];
                    znew = znew - 2.0 * time_to_wall_y * Vz[i];
                    Rz_artificial[i] -= 2.0 * time_to_wall_y * Vz[i];

                    Vx[i] = -Vx[i];
                    Vy[i] = -Vy[i];
                    Vz[i] = -Vz[i];
                }
                // Re-apply periodic boundaries for x and z after potential bounce-back
                if (xnew >= box_x || xnew < 0) xnew = mod(xnew, L_x);
                if (open == 0) { // Periodic boundaries in z
                    if (znew < 0 || znew >= box_z) znew = mod(znew, L_z);
                } else { // Open system in z
                    if (znew < 0) znew = 0.0;
                    if (znew >= box_z) znew = box_z - 0.0001;
                }
                if (ynew >= box_y || ynew < 0) ynew = mod(ynew, L_y); // This line seems redundant if solid walls are handled above
            }
        }

        // Update positions if not adsorbed
        if (catalyst[i][1] == 0 || catalyst[i][1] == 3) {
            Rx[i] = xnew;
            Ry[i] = ynew;
        }
        Rz[i] = znew; // Always update z-position

        // --- CONV_DIFF specific logic (particle re-insertion/marking) ---
        if (CONV_DIFF == 1) { // Initial rectangle concentration method
            if (znew < alpha * box_z) {
                double rand_i = alpha / 2.0 * box_z * (double)rand() / (double)RAND_MAX;
                if (znew < alpha / 2.0 * box_z) {
                    if (znew > rand_i) {
                        double randi = (double)rand() / (double)RAND_MAX;
                        if (randi < (double)x1 * beta / (double)gamma1)
                            Ry[i] = beta * box_y * (double)rand() / (double)RAND_MAX;
                        else if (randi < ((double)x1 * beta / (double)gamma1 + (1.0 - 2.0 * beta) * (double)x2 / (double)gamma1))
                            Ry[i] = beta * box_y + (1.0 - 2.0 * beta) * box_y * (double)rand() / (double)RAND_MAX;
                        else
                            Ry[i] = (1.0 - beta) * box_y + beta * box_y * (double)rand() / (double)RAND_MAX;
                    }
                }
                if (znew > alpha / 2.0 * box_z) {
                    if (znew < alpha * box_z - rand_i) {
                        double randi = (double)rand() / (double)RAND_MAX;
                        if (randi < (double)x1 * beta / (double)gamma1)
                            Ry[i] = beta * box_y * (double)rand() / (double)RAND_MAX;
                        else if (randi < ((double)x1 * beta / (double)gamma1 + (1.0 - 2.0 * beta) * (double)x2 / (double)gamma1))
                            Ry[i] = beta * box_y + (1.0 - 2.0 * beta) * box_y * (double)rand() / (double)RAND_MAX;
                        else
                            Ry[i] = (1.0 - beta) * box_y + beta * box_y * (double)rand() / (double)RAND_MAX;
                    }
                }
                velocity_changer(i); // Update velocity based on flow profile
            }
        } else if (CONV_DIFF == 2 && CONV_DIFF_mass == 0) { // Axial dispersion, same mass particles
            if (t >= tmin && (double)((int)(t / DT) % (int)(t_interval / DT)) < (t_step / DT)) {
                if (znew < alpha * box_z) {
                    double randi = (double)rand() / (double)RAND_MAX;
                    if (randi < 0.5 && activate[i][0] == 0) {
                        activate[i][0] = -1; // Mark for removal/change
                    } else if (randi < 0.5 && activate[i][0] == -1) {
                        activate[i][0] = 0; // Mark as active
                    }
                }
                Rz_artificial[i] = Rz[i]; // Reset artificial z-position
            }
        } else if (CONV_DIFF == 2 && CONV_DIFF_mass == 1) { // Axial dispersion, different mass particles
            if (t >= tmin && (double)((int)(t / DT) % (int)(t_interval / DT)) < (t_step / DT)) {
                if (znew < alpha * box_z) {
                    double randi = (double)rand() / (double)RAND_MAX;
                    if (randi < 0.5 && activate[i][0] == 0) {
                        activate[i][0] = -1; // Mark for removal/change
                        // Add to mass-related linked list
                        int cellx_mass = (int)Rx[i];
                        int celly_mass = (int)Ry[i];
                        int cellz_mass = (int)znew;
                        list_mass[i] = head_mass[cellx_mass][celly_mass][cellz_mass];
                        head_mass[cellx_mass][celly_mass][cellz_mass] = i;
                    } else if (randi < 0.5 && activate[i][0] == -1) {
                        activate[i][0] = 0; // Mark as active
                    }
                }
                Rz_artificial[i] = Rz[i]; // Reset artificial z-position
            }
        }
    }
}

/**
 * @brief Performs the collision step for SRD particles within cells.
 * Includes ghost particles for solid walls and applies a Galilean thermostat.
 */
void collide_particles()
{
    int xcell, ycell, zcell, i, j, k;
    double avg_vel[3], rel_vel[3], rot_vel[3], rot_matrix[3][3];
    double axis_z, axis_theta, axis_phi;
    double cell_mass_sum, cell_particle_count;
    double kinetic_energy_rel = 0.0;
    int num_free_particles = 0;
    double scale_factor = 1.0;

    // Calculate average mass density for ghost particle insertion
    double average_mass_density = ((double)(conc * massA + (1.0 - conc)) * Npart) / (box_x * box_y * box_z);

    if (therm == 1) { // If thermostat is active, calculate scaling factor
        scale_factor = sqrt(1.0 + 1.0 / 2.0 * ((kTo / kTmeas) - 1.0));
    }

    for (xcell = 0; xcell < L_x; xcell++) {
        for (ycell = 0; ycell < L_y + 1; ycell++) { // L_y + 1 to include wall cells
            for (zcell = 0; zcell < L_z; zcell++) {
                avg_vel[0] = 0.0; // Reset average velocity for current cell
                avg_vel[1] = 0.0;
                avg_vel[2] = 0.0;
                cell_mass_sum = 0.0;
                cell_particle_count = 0.0;

                i = head[xcell][ycell][zcell];
                while (i >= 0) {
                    if (i >= Nsrd) { // Foreign particles
                        cell_mass_sum += massA;
                        avg_vel[0] += massA * Vx[i];
                        avg_vel[1] += massA * Vy[i];
                        avg_vel[2] += massA * Vz[i];
                    } else { // SRD particles
                        cell_mass_sum += 1.0;
                        avg_vel[0] += Vx[i];
                        avg_vel[1] += Vy[i];
                        avg_vel[2] += Vz[i];
                    }
                    cell_particle_count += 1.0;
                    i = list[i];
                }

                // Ghost particle insertion for solid walls (Lamura & Gompper)
                if (walls == 1 && nearwalls[xcell][ycell][zcell] == 1) {
                    double vel_multiplier = sqrt(kTo);
                    // Add ghost particles if cell mass is less than average
                    while (cell_mass_sum < average_mass_density) {
                        avg_vel[0] += vel_multiplier * sqrt(average_mass_density - cell_mass_sum) * neargauss();
                        avg_vel[1] += vel_multiplier * sqrt(average_mass_density - cell_mass_sum) * neargauss();
                        avg_vel[2] += vel_multiplier * sqrt(average_mass_density - cell_mass_sum) * neargauss();
                        cell_mass_sum = average_mass_density; // Cell mass is now effectively average
                    }
                    // Remove ghost particles if cell mass is more than average (original logic had this, but seems to be a bug)
                    // The original code had a similar loop for mass > averageM, which would subtract velocities.
                    // This part is usually for momentum conservation, not mass. Keeping it commented for now.
                    /*
                    while (cell_mass_sum > average_mass_density) {
                        avg_vel[0] -= vel_multiplier * sqrt(cell_mass_sum - average_mass_density) * neargauss();
                        avg_vel[1] -= vel_multiplier * sqrt(cell_mass_sum - average_mass_density) * neargauss();
                        avg_vel[2] -= vel_multiplier * sqrt(cell_mass_sum - average_mass_density) * neargauss();
                        cell_mass_sum = average_mass_density;
                    }
                    */
                }

                if (cell_particle_count > 1.0 || (walls == 1 && nearwalls[xcell][ycell][zcell] == 1 && cell_particle_count > 0)) {
                    // Calculate average velocity for the cell
                    avg_vel[0] /= cell_mass_sum;
                    avg_vel[1] /= cell_mass_sum;
                    avg_vel[2] /= cell_mass_sum;

                    // Choose random rotation axis (ax,ay,az)
                    axis_z = 2.0 * (double)rand() / (double)RAND_MAX - 1.0;
                    axis_theta = sqrt(1.0 - axis_z * axis_z);
                    axis_phi = 2.0 * M_PI * (double)rand() / (double)RAND_MAX;
                    double axis_x = axis_theta * cos(axis_phi);
                    double axis_y = axis_theta * sin(axis_phi);

                    // Rotation matrix (fixed 90 degrees around axis ax,ay,az)
                    // This is a rotation by angle M_PI/2 (90 degrees)
                    rot_matrix[0][0] = axis_x * axis_x;
                    rot_matrix[0][1] = axis_x * axis_y - axis_z;
                    rot_matrix[0][2] = axis_x * axis_z + axis_y;
                    rot_matrix[1][0] = axis_x * axis_y + axis_z;
                    rot_matrix[1][1] = axis_y * axis_y;
                    rot_matrix[1][2] = axis_y * axis_z - axis_x;
                    rot_matrix[2][0] = axis_x * axis_z - axis_y;
                    rot_matrix[2][1] = axis_y * axis_z + axis_x;
                    rot_matrix[2][2] = axis_z * axis_z;

                    num_free_particles--; // Subtract 1 particle because c.o.m. velocity is fixed

                    i = head[xcell][ycell][zcell];
                    while (i >= 0) { // Rotate each velocity relative to average
                        rel_vel[0] = Vx[i] - avg_vel[0];
                        rel_vel[1] = Vy[i] - avg_vel[1];
                        rel_vel[2] = Vz[i] - avg_vel[2];

                        if (therm == 1) { // Apply thermostat on relative velocities
                            rel_vel[0] *= scale_factor;
                            rel_vel[1] *= scale_factor;
                            rel_vel[2] *= scale_factor;

                            // Calculate new global temperature
                            if (i >= Nsrd) { // Foreign particles
                                kinetic_energy_rel += massA * (rel_vel[0] * rel_vel[0] + rel_vel[1] * rel_vel[1] + rel_vel[2] * rel_vel[2]);
                            } else { // SRD particles
                                kinetic_energy_rel += rel_vel[0] * rel_vel[0] + rel_vel[1] * rel_vel[1] + rel_vel[2] * rel_vel[2];
                            }
                            num_free_particles++;
                        }

                        // Apply rotation matrix
                        for (j = 0; j < 3; j++) {
                            rot_vel[j] = 0;
                            for (k = 0; k < 3; k++) {
                                rot_vel[j] += rot_matrix[j][k] * rel_vel[k];
                            }
                        }

                        Vx[i] = avg_vel[0] + rot_vel[0];
                        Vy[i] = avg_vel[1] + rot_vel[1];
                        Vz[i] = avg_vel[2] + rot_vel[2];

                        i = list[i];
                    } // End while (i>=0)
                } // End if (cell_particle_count > 1.0 ...)
            }
        }
    } // End loop over cells

    if (therm == 1 && num_free_particles > 0) {
        kTmeas = kinetic_energy_rel / (3.0 * num_free_particles); // Update measured temperature
    }
}

/**
 * @brief Generates two Gaussian distributed random numbers using the Box-Muller transform.
 * @param gauss1 Pointer to store the first Gaussian random number.
 * @param gauss2 Pointer to store the second Gaussian random number.
 */
void GAUSS(double *gauss1, double *gauss2)
{
    double x1, x2, twou, radius, theta;

    x1 = (double)rand() / (double)RAND_MAX;
    x2 = (double)rand() / (double)RAND_MAX;
    twou = -2 * log(1.0 - x1);
    radius = sqrt(twou);
    theta = 2 * M_PI * x2;
    *gauss1 = radius * cos(theta);
    *gauss2 = radius * sin(theta);
}

/**
 * @brief Generates a nearly Gaussian distributed random number of unit variance
 * and average of 0 by summing 12 uniform random numbers (Central Limit Theorem).
 * @return A single Gaussian-like random number.
 */
double neargauss()
{
    int i;
    double sum = -6.0; // Sum of 12 uniform random numbers (0-1) has mean 6.0, so subtract 6.0 for mean 0.
    for (i = 0; i < 12; i++) {
        sum += (double)rand() / (double)RAND_MAX;
    }
    return sum;
}

/**
 * @brief Placeholder function for wall velocity.
 * Original code had this, but its usage is unclear without more context.
 * @return A random value based on a Gaussian distribution.
 */
double wall_vel()
{
    double random = (double)rand() / (double)RAND_MAX;
    double X = sqrt(-2 * log(1 - random));
    return X;
}

/**
 * @brief Returns the nearest integer for a double number.
 * Handles both positive and negative numbers.
 * @param r The double number.
 * @return The nearest integer.
 */
int nint(double r)
{
    return (r > 0.0) ? (int)(r + 0.5) : (int)(r - 0.5);
}

/**
 * @brief Adjusts the center of mass velocity to zero.
 * The original implementation was commented out or had complex logic.
 * This version is a simplified placeholder based on the original comments.
 * It iterates through particles and adjusts their velocities to ensure
 * the overall center of mass velocity is zero.
 */
void COM_vel_zero()
{
    double sum_vx = 0.0;
    double sum_vy = 0.0;
    double sum_vz = 0.0;
    double total_mass = 0.0;

    for (int i = 0; i < Nsrd; i++) { // SRD particles (mass 1.0)
        sum_vx += Vx[i];
        sum_vy += Vy[i];
        sum_vz += Vz[i];
        total_mass += 1.0;
    }
    for (int i = Nsrd; i < Npart; i++) { // Foreign particles (mass massA)
        sum_vx += massA * Vx[i];
        sum_vy += massA * Vy[i];
        sum_vz += massA * Vz[i];
        total_mass += massA;
    }

    if (total_mass > 0) {
        sum_vx /= total_mass;
        sum_vy /= total_mass;
        sum_vz /= total_mass;

        for (int i = 0; i < Npart; i++) {
            if (i < Nsrd) {
                Vx[i] -= sum_vx;
                Vy[i] -= sum_vy;
                Vz[i] -= sum_vz;
            } else {
                Vx[i] -= sum_vx / massA;
                Vy[i] -= sum_vy / massA;
                Vz[i] -= sum_vz / massA;
            }
        }
    }
}

/**
 * @brief Handles adsorption, desorption, and reaction processes on the catalytic surface.
 * This function updates particle states (adsorbed/bulk) and surface coverages.
 */
void adsorption()
{
    // Debug print (original code had this, keeping it for now)
    if (mod(t, 100) < 0.2) {
        printf("%d\t%d\t%f\t%f\t%f\t%d\n", adsA, adsB, thetaA, thetaB, theta, adsA + adsB + NA + Nsrd);
    }

    // --- Desorption and Reaction of Adsorbed A Particles ---
    for (int i = 0; i < adsA; i++) {
        double random_val = (double)rand() / (double)RAND_MAX;
        if (random_val < p_desA) { // Desorption of A
            thetaA -= 1.0 / (double)catalytic_sites;
            theta = 1 - thetaA - thetaB;
            Nsrd++; // Adsorbed A becomes an SRD particle
            Npart++;

            // Shift particles to make space for the new SRD particle
            for (int k = Npart - 1; k >= Nsrd; k--) {
                Rx[k] = Rx[k - 1];
                Ry[k] = Ry[k - 1];
                Rz[k] = Rz[k - 1];
                Vx[k] = Vx[k - 1];
                Vy[k] = Vy[k - 1];
                Vz[k] = Vz[k - 1];
                catalyst[k][0] = catalyst[k - 1][0];
                catalyst[k][1] = catalyst[k - 1][1];
            }

            // Assign new velocities and position to the desorbed particle
            Vx[Nsrd - 1] = neargauss();
            Vy[Nsrd - 1] = neargauss();
            if (Vy[Nsrd - 1] < 0.0) Vy[Nsrd - 1] = -Vy[Nsrd - 1]; // Ensure positive y-velocity for desorption
            Vz[Nsrd - 1] = neargauss();

            // Position the desorbed particle near the surface
            Rx[Nsrd - 1] = mod(Rx_adsA[i] + Vx[Nsrd - 1] * random_val * DT, L_x);
            Ry[Nsrd - 1] = mod(Ry_adsA[i] + Vy[Nsrd - 1] * random_val * DT, L_y);
            Rz[Nsrd - 1] = mod(Rz_adsA[i] + Vz[Nsrd - 1] * random_val * DT, L_z);

            catalyst[Nsrd - 1][0] = 0; // Mark as SRD particle
            catalyst[Nsrd - 1][1] = 0; // Mark as not adsorbed

            // Remove the adsorbed particle from the adsA list
            for (int k = i; k < adsA - 1; k++) { // Bug fix: k++ instead of k--
                Rx_adsA[k] = Rx_adsA[k + 1];
                Ry_adsA[k] = Ry_adsA[k + 1];
                Rz_adsA[k] = Rz_adsA[k + 1];
            }
            adsA--;
            i--; // Adjust loop index as an element was removed
        } else if (random_val < p_desA + p_react) { // Reaction of A to B
            thetaB += 1.0 / (double)catalytic_sites;
            thetaA -= 1.0 / (double)catalytic_sites;

            // Move particle from adsA list to adsB list
            Rx_adsB[adsB] = Rx_adsA[i];
            Ry_adsB[adsB] = Ry_adsA[i];
            Rz_adsB[adsB++] = Rz_adsA[i];

            // Remove from adsA list
            for (int k = i; k < adsA - 1; k++) { // Bug fix: k++ instead of k--
                Rx_adsA[k] = Rx_adsA[k + 1];
                Ry_adsA[k] = Ry_adsA[k + 1];
                Rz_adsA[k] = Rz_adsA[k + 1];
            }
            adsA--;
            i--; // Adjust loop index
            theta = 1 - thetaA - thetaB;
        }
    }

    // --- Desorption of Adsorbed B Particles ---
    for (int i = 0; i < adsB; i++) {
        double random_val = (double)rand() / (double)RAND_MAX;
        if (random_val < p_desB) { // Desorption of B
            thetaB -= 1.0 / (double)catalytic_sites;
            theta = 1 - thetaA - thetaB;
            NA++; // Adsorbed B becomes a foreign particle
            Npart++;

            // Assign new velocities and position to the desorbed particle
            Vx[Npart - 1] = neargauss();
            Vy[Npart - 1] = neargauss();
            if (Vy[Npart - 1] < 0.0) Vy[Npart - 1] = -Vy[Npart - 1];
            Vz[Npart - 1] = neargauss();

            // Position the desorbed particle near the surface
            Rx[Npart - 1] = mod(Rx_adsB[i] + Vx[Npart - 1] * random_val * DT, L_x);
            Ry[Npart - 1] = mod(Ry_adsB[i] + Vy[Npart - 1] * random_val * DT, L_y);
            Rz[Npart - 1] = mod(Rz_adsB[i] + Vz[Npart - 1] * random_val * DT, L_z);

            catalyst[Npart - 1][0] = 3; // Mark as foreign particle
            catalyst[Npart - 1][1] = 3; // Mark as not adsorbed

            // Remove the adsorbed particle from the adsB list
            for (int k = i; k < adsB - 1; k++) { // Bug fix: k++ instead of k--
                Rx_adsB[k] = Rx_adsB[k + 1];
                Ry_adsB[k] = Ry_adsB[k + 1];
                Rz_adsB[k] = Rz_adsB[k + 1];
            }
            adsB--;
            i--; // Adjust loop index
        }
    }

    // --- Update Catalyst States for Particles in Bulk (if they were marked as adsorbed) ---
    // This section seems to handle particles that were marked as adsorbed but are now in bulk
    // due to previous desorption/reaction events.
    for (int i = 0; i < Npart; i++) {
        if (catalyst[i][0] == 0 && catalyst[i][1] == 1) { // SRD particle, marked as adsorbed A
            // This particle should be moved to adsA list
            Rx_adsA[adsA] = Rx[i];
            Ry_adsA[adsA] = Ry[i];
            Rz_adsA[adsA++] = Rz[i];

            // Remove from bulk particle list
            for (int k = i; k < Npart - 1; k++) { // Bug fix: k++ instead of k--
                Rx[k] = Rx[k + 1]; Vx[k] = Vx[k + 1];
                Ry[k] = Ry[k + 1]; Vy[k] = Vy[k + 1];
                Rz[k] = Rz[k + 1]; Vz[k] = Vz[k + 1];
                catalyst[k][0] = catalyst[k + 1][0];
                catalyst[k][1] = catalyst[k + 1][1];
            }
            Nsrd--;
            Npart--;
            thetaA += 1.0 / (double)catalytic_sites;
            theta = 1 - thetaA - thetaB;
            i--; // Adjust loop index
            continue;
        }
        if (catalyst[i][0] == 0 && catalyst[i][1] == 2) { // SRD particle, marked as adsorbed B (after reaction)
            // This particle should be moved to adsB list
            Rx_adsB[adsB] = Rx[i];
            Ry_adsB[adsB] = Ry[i];
            Rz_adsB[adsB++] = Rz[i];

            // Remove from bulk particle list
            for (int k = i; k < Npart - 1; k++) { // Bug fix: k++ instead of k--
                Rx[k] = Rx[k + 1]; Vx[k] = Vx[k + 1];
                Ry[k] = Ry[k + 1]; Vy[k] = Vy[k + 1];
                Rz[k] = Rz[k + 1]; Vz[k] = Vz[k + 1];
                catalyst[k][0] = catalyst[k + 1][0];
                catalyst[k][1] = catalyst[k + 1][1];
            }
            Nsrd--;
            Npart--;
            thetaB += 1.0 / (double)catalytic_sites;
            theta = 1 - thetaA - thetaB;
            i--; // Adjust loop index
            continue;
        }
        if (catalyst[i][0] == 0 && catalyst[i][1] == 3) { // SRD particle, marked as desorbed B (after reaction and desorption)
            // This particle should become a foreign particle in bulk
            double Vxtemp = -Vx[i];
            double Vytemp = -Vy[i];
            double Vztemp = -Vz[i];
            double Rxtemp = Rx[i] + Vxtemp * timey; // 'timey' is a global variable, its usage here is unclear
            double Rytemp = Ry[i] + Vytemp * timey;
            double Rztemp = Rz[i] + Vztemp * timey;

            // Remove from current position and re-add as foreign particle at the end
            for (int k = i; k < Npart - 1; k++) { // Bug fix: k++ instead of k--
                Rx[k] = Rx[k + 1]; Vx[k] = Vx[k + 1];
                Ry[k] = Ry[k + 1]; Vy[k] = Vy[k + 1];
                Rz[k] = Rz[k + 1]; Vz[k] = Vz[k + 1];
                catalyst[k][0] = catalyst[k + 1][0];
                catalyst[k][1] = catalyst[k + 1][1];
            }
            Nsrd--;
            NA++;
            Rx[Npart - 1] = Rxtemp;
            Ry[Npart - 1] = Rytemp;
            Rz[Npart - 1] = Rztemp;
            Vx[Npart - 1] = Vxtemp;
            Vy[Npart - 1] = Vytemp;
            Vz[Npart - 1] = Vztemp;
            catalyst[Npart - 1][0] = 3; // Mark as foreign particle
            catalyst[Npart - 1][1] = 3; // Mark as not adsorbed
            i--; // Adjust loop index
            continue;
        }
        if (catalyst[i][0] == 3 && catalyst[i][1] == 2) { // Foreign particle, marked as adsorbed B
            // This particle should be moved to adsB list
            Rx_adsB[adsB] = Rx[i];
            Ry_adsB[adsB] = Ry[i];
            Rz_adsB[adsB++] = Rz[i];

            // Remove from bulk particle list
            for (int k = i; k < Npart - 1; k++) { // Bug fix: k++ instead of k--
                Rx[k] = Rx[k + 1]; Vx[k] = Vx[k + 1];
                Ry[k] = Ry[k + 1]; Vy[k] = Vy[k + 1];
                Rz[k] = Rz[k + 1]; Vz[k] = Vz[k + 1];
                catalyst[k][0] = catalyst[k + 1][0];
                catalyst[k][1] = catalyst[k + 1][1];
            }
            NA--;
            Npart--;
            thetaB += 1.0 / (double)catalytic_sites;
            theta = 1 - thetaA - thetaB;
            i--; // Adjust loop index
            continue;
        }
    }
}