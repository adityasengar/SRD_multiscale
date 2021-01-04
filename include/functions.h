#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

/* main.c */
void initial_conditions();
// void processInput(int argn, char * argv[], int stepmd, int print_energy, int trajoutput); // Not used in main.c, removed

/* io.c */
void record_trajectories();
void energy_momentum( int info, double time );
void updatecorrelator_MSD(int reset, double massA);
void velocityprofile(int reset, double time);
void slit_tracker(int reset, double time);
void profile_maker(int reset, double time);
void velocity_correlator(int reset, double dt);
void dual_velocity_correlator(int reset, double massA);
void new_diff(int reset, double t);
void something(int reset, double massA); // Renamed from original 'something' for clarity if possible, but keeping for now.

/* calc.c */
double mod(double r, int a);
void createlinkedlist();
// void calculateforces(); // Commented out in original calc.c, so not including prototype
void updatevelocities();
double min3(double r1,double r2,double r3);
void updatepositions();
void collide_particles();
double neargauss();
int nint( double r);
void GAUSS( double *gauss1, double *gauss2);
void COM_vel_zero(); // Added prototype for COM_vel_zero
void adsorption(); // Moved from velocity.c to calc.c

/* boundaries.c */
void determinenearwalls();

/* memory.c */
void allocateMemory();
void deallocateMemory();
void alloc_memory_ACF(int Npart, int numcor);
void deallocateMemory_ACF();
void alloc_memory_CONV_DIFF_mass(int Npart);

/* velocity.c */
void homogenous(int reset, double tx);
void update_velz_bin(int reset, int total_time);
void velocity_changer(int i);

#endif /* FUNCTIONS_H_ */
