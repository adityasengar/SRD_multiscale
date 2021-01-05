#include <stdio.h>
#include <stdlib.h>
#include "variables.h" // Assuming variables.h defines global arrays like Rx, Ry, etc.
#include "constants.h" // Assuming constants.h defines Npartmax, Lboxmax etc.
#include "memory.h"

// Function to allocate memory for all global arrays
void allocateMemory() {
    // Allocate memory for SRD particle positions and velocities
    Rx = (double*)malloc(Npart * sizeof(double));
    Ry = (double*)malloc(Npart * sizeof(double));
    Rz = (double*)malloc(Npart * sizeof(double));
    Vx = (double*)malloc(Npart * sizeof(double));
    Vy = (double*)malloc(Npart * sizeof(double));
    Vz = (double*)malloc(Npart * sizeof(double));
    Rz_artificial = (double*)malloc(Npart * sizeof(double));

    // Allocate memory for linked-list
    list = (int*)malloc(Npart * sizeof(int));

    // Allocate memory for head and nearwalls 3D arrays
    head = (int***)malloc(L_x * sizeof(int**));
    nearwalls = (int***)malloc(L_x * sizeof(int**));
    catalyst = (int**)malloc(Npart * sizeof(int*));
    Rx_adsA = (double*)malloc(Npart * sizeof(double));
    Ry_adsA = (double*)malloc(Npart * sizeof(double));
    Rz_adsA = (double*)malloc(Npart * sizeof(double));
    Rx_adsB = (double*)malloc(Npart * sizeof(double));
    Ry_adsB = (double*)malloc(Npart * sizeof(double));
    Rz_adsB = (double*)malloc(Npart * sizeof(double));


    for (int i = 0; i < L_x; i++) {
        head[i] = (int**)malloc((L_y + 1) * sizeof(int*));
        nearwalls[i] = (int**)malloc((L_y + 1) * sizeof(int*));
        for (int j = 0; j < (L_y + 1); j++) {
            head[i][j] = (int*)malloc((L_z + 1) * sizeof(int)); // L_z + 1 for open boundary
            nearwalls[i][j] = (int*)malloc(L_z * sizeof(int));
        }
    }

    for (int i = 0; i < Npart; i++) {
        catalyst[i] = (int*)malloc(2 * sizeof(int));
    }

    // Initialize arrays (optional, but good practice)
    for (int i = 0; i < Npart; i++) {
        Rx[i] = Ry[i] = Rz[i] = 0.0;
        Vx[i] = Vy[i] = Vz[i] = 0.0;
        Rz_artificial[i] = 0.0;
        list[i] = -1;
        catalyst[i][0] = 0; // Default to not a catalyst particle
        catalyst[i][1] = 0; // Default history
    }

    for (int i = 0; i < L_x; i++) {
        for (int j = 0; j < (L_y + 1); j++) {
            for (int k = 0; k < (L_z + 1); k++) { // L_z + 1 for open boundary
                head[i][j][k] = -1;
            }
            for (int k = 0; k < L_z; k++) {
                nearwalls[i][j][k] = 0;
            }
        }
    }

    // Allocate memory for velocity profile arrays if needed
    vel_sum = (double*)calloc(ybin, sizeof(double));
    counter = (int*)calloc(ybin, sizeof(int));
    y_sum = (double*)calloc(ybin, sizeof(double));
    z_sum = (double*)calloc(zbin, sizeof(double));
    z_sum1 = (double*)calloc(zbin, sizeof(double));
    velz_meas = (double*)calloc(velz_bin, sizeof(double));
    velz_flow = (double*)calloc(velz_bin, sizeof(double));
    N_average = (double*)calloc(velz_bin, sizeof(double));
}

// Function to deallocate memory for all global arrays
void deallocateMemory() {
    free(Rx);
    free(Ry);
    free(Rz);
    free(Vx);
    free(Vy);
    free(Vz);
    free(Rz_artificial);
    free(list);

    for (int i = 0; i < L_x; i++) {
        for (int j = 0; j < (L_y + 1); j++) {
            free(head[i][j]);
            free(nearwalls[i][j]);
        }
        free(head[i]);
        free(nearwalls[i]);
    }
    free(head);
    free(nearwalls);

    for (int i = 0; i < Npart; i++) {
        free(catalyst[i]);
    }
    free(catalyst);
    free(Rx_adsA);
    free(Ry_adsA);
    free(Rz_adsA);
    free(Rx_adsB);
    free(Ry_adsB);
    free(Rz_adsB);

    free(vel_sum);
    free(counter);
    free(y_sum);
    free(z_sum);
    free(z_sum1);
    free(velz_meas);
    free(velz_flow);
    free(N_average);
}

// Placeholder for specific allocation functions
int** int_2D_array(int dim1, int dim2) {
    int** arr = (int**)malloc(dim1 * sizeof(int*));
    for (int i = 0; i < dim1; i++) {
        arr[i] = (int*)malloc(dim2 * sizeof(int));
    }
    return arr;
}

int* int_1D_array(int dim1) {
    return (int*)malloc(dim1 * sizeof(int));
}

int*** int_3D_array(int dim1, int dim2, int dim3) {
    int*** arr = (int***)malloc(dim1 * sizeof(int**));
    for (int i = 0; i < dim1; i++) {
        arr[i] = (int**)malloc(dim2 * sizeof(int*));
        for (int j = 0; j < dim2; j++) {
            arr[i][j] = (int*)malloc(dim3 * sizeof(int));
        }
    }
    return arr;
}

void free_1D_array(void* ptr) {
    free(ptr);
}

void free_2D_array(void** ptr, int dim1) {
    for (int i = 0; i < dim1; i++) {
        free(ptr[i]);
    }
    free(ptr);
}

void free_3D_array(void*** ptr, int dim1, int dim2) {
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            free(ptr[i][j]);
        }
        free(ptr[i]);
    }
    free(ptr);
}

void alloc_memory_ACF(int Npart_val, int numcor_val) {
    cor = (double*)calloc(numcor_val, sizeof(double));
    cor1 = (double*)calloc(numcor_val, sizeof(double));
    cor2 = (double*)calloc(numcor_val, sizeof(double));
    ACFcount = (int*)calloc(numcor_val, sizeof(int));
    ACFcount1 = (int*)calloc(numcor_val, sizeof(int));
    ACFcount2 = (int*)calloc(numcor_val, sizeof(int));

    storex = (double**)malloc(Npart_val * sizeof(double*));
    storey = (double**)malloc(Npart_val * sizeof(double*));
    storez = (double**)malloc(Npart_val * sizeof(double*));
    storex1 = (double**)malloc(Npart_val * sizeof(double*));
    storey1 = (double**)malloc(Npart_val * sizeof(double*));
    storez1 = (double**)malloc(Npart_val * sizeof(double*));
    storex2 = (double**)malloc(Npart_val * sizeof(double*));
    storey2 = (double**)malloc(Npart_val * sizeof(double*));
    storez2 = (double**)malloc(Npart_val * sizeof(double*));

    for (int i = 0; i < Npart_val; i++) {
        storex[i] = (double*)malloc(numcor_val * sizeof(double));
        storey[i] = (double*)malloc(numcor_val * sizeof(double));
        storez[i] = (double*)malloc(numcor_val * sizeof(double));
        storex1[i] = (double*)malloc(numcor_val * sizeof(double));
        storey1[i] = (double*)malloc(numcor_val * sizeof(double));
        storez1[i] = (double*)malloc(numcor_val * sizeof(double));
        storex2[i] = (double*)malloc(numcor_val * sizeof(double));
        storey2[i] = (double*)malloc(numcor_val * sizeof(double));
        storez2[i] = (double*)malloc(numcor_val * sizeof(double));
    }
}

void deallocateMemory_ACF() {
    free(cor);
    free(cor1);
    free(cor2);
    free(ACFcount);
    free(ACFcount1);
    free(ACFcount2);

    for (int i = 0; i < Npart; i++) {
        free(storex[i]);
        free(storey[i]);
        free(storez[i]);
        free(storex1[i]);
        free(storey1[i]);
        free(storez1[i]);
        free(storex2[i]);
        free(storey2[i]);
        free(storez2[i]);
    }
    free(storex);
    free(storey);
    free(storez);
    free(storex1);
    free(storey1);
    free(storez1);
    free(storex2);
    free(storey2);
    free(storez2);
}

void alloc_memory_CONV_DIFF_mass(int Npart_val) {
    list_mass = int_1D_array(Npart_val);
    head_mass = int_3D_array(LBOX_MAX, LBOX_MAX + 1, LBOX_MAX); // Assuming Lboxmax is defined in constants.h
}
