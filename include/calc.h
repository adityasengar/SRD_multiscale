#ifndef CALC_H_INCLUDED
#define CALC_H_INCLUDED

// Function prototypes for calc.c
double mod(double r, int a);
void createlinkedlist();
void updatevelocities();
double min3(double r1,double r2,double r3);
void updatepositions();
void collide_particles();
double neargauss();
int nint( double r);
void GAUSS( double *gauss1, double *gauss2);
void COM_vel_zero(); // Added prototype for COM_vel_zero
void adsorption(); // Added prototype for adsorption

#endif // CALC_H_INCLUDED
