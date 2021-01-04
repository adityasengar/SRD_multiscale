#ifndef MEMORY_H_
#define MEMORY_H_

// Function declarations for memory management
void allocateMemory();
void deallocateMemory();

// Placeholder for specific allocation functions if they exist
// These are inferred from usage in calc.c and main.c
int** int_2D_array(int dim1, int dim2);
int* int_1D_array(int dim1);
int*** int_3D_array(int dim1, int dim2, int dim3);
void free_1D_array(void* ptr);
void free_2D_array(void** ptr, int dim1);
void free_3D_array(void*** ptr, int dim1, int dim2);
void alloc_memory_ACF(int Npart, int numcor);
void deallocateMemory_ACF();
void alloc_memory_CONV_DIFF_mass(int Npart);

#endif /* MEMORY_H_ */
