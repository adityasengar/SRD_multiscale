#ifndef VELOCITY_H_
#define VELOCITY_H_

#include <stdio.h>
#include "constants.h"
#include "variables.h"

// Function prototypes for velocity.c
void homogenous(int reset, double tx);
void update_velz_bin(int reset, int total_time);
void velocity_changer(int i);

#endif /* VELOCITY_H_ */
